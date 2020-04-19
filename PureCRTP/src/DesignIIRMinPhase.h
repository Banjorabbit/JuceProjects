#pragma once
#include "BaseClasses/PureCRTP.h"
#include "FFT.h"
#include "FrequencyDomain/MinPhaseSpectrum.h"
#include <unsupported/Eigen/Polynomials>

struct O::DesignIIRMinPhase
{
	Real2D SOS;
	Float Gain;
};

class DesignIIRMinPhase : public Base<DesignIIRMinPhase, I::Real2D, O::DesignIIRMinPhase>
{
	friend Base<DesignIIRMinPhase, I::Real2D, O::DesignIIRMinPhase>;

public:
	static int GetNSOS(int N) { return static_cast<int>(std::ceil(static_cast<float>(N) / 2)); } // get number of SOS

private:
	FFTReal FFT;
	MinPhaseSpectrum MinPhaseCalculation;

	struct Coefficients
	{
		int NBands = 65; // NBands should be significantly higher than C.N
		float SampleRate = 44100.f;
		ConstrainedType<int> N = { 8, 1, 32 }; // very high orders are not numerically robust when calculating zeros/poles (roots)
		enum WeightTypes { Linear, MelScale};
		WeightTypes WeightType = MelScale;
	} C;

	struct Parameters { } P;

	struct Data
	{
		int NA, NB;
		Eigen::PolynomialSolver<float, Eigen::Dynamic> Polynomialsolver;
		int FFTSize;
		Eigen::MatrixXf R;
		Eigen::ArrayXf Weight;
		void Reset() {  }
		bool InitializeMemory(const Coefficients& c)
		{
			NB = c.N;
			NA = c.N;
			FFTSize = (c.NBands - 1) * 2;
			R.resize(NA + NB + 1, NA + NB + 1);
			R.block(NA, NA, NB + 1, NB + 1).setIdentity();
			Weight.resize(c.NBands);
			switch (c.WeightType)
			{
			case Coefficients::Linear:
				Weight.setOnes();
				break;
			case Coefficients::MelScale:
				Weight(0) = 1.6f;
				Eigen::ArrayXf freq = Eigen::ArrayXf::LinSpaced(c.NBands - 1, 1, static_cast<float>(c.NBands - 1))*c.SampleRate/FFTSize;
				Weight.tail(c.NBands - 1) = 2595.f*(1.f + freq / 700).log10() / freq;
				break;
			}
			return true;
		}
		size_t GetAllocatedMemorySize() const
		{
			auto size = R.GetAllocatedMemorySize();
			size += Weight.GetAllocatedMemorySize();
			return size;
		}
		void OnParameterChange(const Parameters& p, const Coefficients& c)	{	}
	} D;

	DEFINEMEMBERALGORITHMS(2, FFT, MinPhaseCalculation)

	auto InitializeMembers()
	{
		auto sFFT = FFT.GetSetup();
		sFFT.Coefficients.FFTSize = D.FFTSize;
		auto flag = FFT.Initialize(sFFT);
		// min phase calculation
		auto setupMPS = MinPhaseCalculation.GetSetup();
		setupMPS.Coefficients.NBands = C.NBands;
		flag &= MinPhaseCalculation.Initialize(setupMPS);
		return flag;
	}

	// convert transfer function to second order sections
	void TF2SOS(const Eigen::Ref<const Eigen::ArrayXf>& tf, Eigen::Ref<Eigen::ArrayXXf> sos, float& gain)
	{
		const auto N = tf.size() - 1;
		const auto N1 = N - 1;

		//// This code is equivalent to the PolynomialSolver (with different ordering of roots), but is slightly slower when profiling. If this is used then Eigen::EigenSolver<Eigen::MatrixXf> needs to be a member in D
		//Eigen::MatrixXf companionMatrix = Eigen::MatrixXf::Zero(N, N);
		//companionMatrix.diagonal(-1) = Eigen::VectorXf::Ones(N - 1);
		//companionMatrix.row(0) = -tf.tail(N).transpose() / tf(0);
		//D.EigenSolver.compute(companionMatrix);
		//Eigen::VectorXcf poles = D.EigenSolver.eigenvalues();

		D.Polynomialsolver.compute(tf.colwise().reverse().matrix());
		Eigen::VectorXcf poles = D.Polynomialsolver.roots();

		Eigen::VectorXcf polesStable(N);
		auto start = 0;
		auto end = N1;
		for (auto i = 0; i < N;i++)
		{
			float magnitudex2 = std::real(poles(i) * std::conj(poles(i)));
			if (poles(i).imag() == 0)
			{
				if (magnitudex2 > 1.f) { polesStable(end) = std::real(poles(i)) / magnitudex2; }
				else { polesStable(end) = poles(i); }
				end--;
			}
			else
			{
				if (magnitudex2 > 1.f) { polesStable(start) = poles(i) / magnitudex2; }
				else { polesStable(start) = poles(i); }
				start++;
			}
		}

		auto nSections = sos.rows();
		sos.col(0).setOnes();
		for (auto i = 0; i < N - 1; i += 2)
		{
			sos(i / 2, 1) = -std::real(polesStable(i) + polesStable(i + 1));
			sos(i / 2, 2) = std::real(polesStable(i) * polesStable(i + 1));
		}
		if (N & 1) // if odd number of roots
		{
			sos(nSections - 1, 1) = -std::real(polesStable(N1));
			sos(nSections - 1, 2) = 0.f;
		}
		float gainSOS = sos.rowwise().sum().prod();
		float gainA = tf.sum();
		gain = gainA / gainSOS;
	}

	void ProcessOn(Input xMag, Output yTime)
	{
		for (auto channel = 0; channel < xMag.cols(); channel++)
		{
			Eigen::ArrayXcf xFreq(C.NBands);
			MinPhaseCalculation.Process(xMag.col(channel), xFreq);

			Eigen::ArrayXf xTime(D.FFTSize);
			Eigen::VectorXf Vd(D.NB + D.NA + 1);

			// lower left and upper right corner
			Eigen::ArrayXcf xPow = xFreq * D.Weight;
			FFT.Inverse(xPow, xTime);
			for (auto i = 0;i < D.NA;i++) // lower left corner
			{
				const auto r2 = std::max(D.NB - i, 0);
				const auto r1 = D.NB + 1 - r2;
				D.R.block(D.NA, i, r1, 1) = -xTime.segment(D.FFTSize - 1 - i, r1);
				D.R.block(D.NA + r1, i, r2, 1) = -xTime.head(r2);
			}
			// upper right corner
			D.R.block(0, D.NA, D.NA, D.NB + 1) = D.R.block(D.NA, 0, D.NB + 1, D.NA).transpose();
			Vd.tail(D.NB + 1) = xTime.head(D.NB + 1);

			// upper left corner
			xPow *= xFreq.conjugate();
			FFT.Inverse(xPow, xTime);
			for (auto i = 0; i < D.NA; i++)
			{
				const auto r2 = D.NA - i;
				const auto r1 = D.NA - r2;
				D.R.block(i, i, r2, 1) = xTime.head(r2);
				D.R.block(0, i, r1, 1) = xTime.segment(1, r1).reverse();
			}
			Vd.head(D.NA) = -xTime.segment(1, D.NA);

			const Eigen::VectorXf th = D.R.llt().solve(Vd);

			Eigen::ArrayXf A(D.NA + 1), B(D.NB+1);
			A(0, channel) = 1.f;
			A.tail(D.NA) = th.head(D.NA);
			B = th.tail(D.NB + 1);

			// convert to sos matrix and gain
			TF2SOS(B, yTime.SOS.leftCols(3), yTime.Gain);
			float gainDen;
			TF2SOS(A, yTime.SOS.rightCols(3), gainDen);
			yTime.Gain /= gainDen;
		}
	}

	void ProcessOff(Input xMag, Output yTime) { yTime.SOS.setZero(); yTime.SOS.col(3).setOnes(); yTime.Gain = 1.f; }
};