#pragma once
#include "../BaseClasses/PureCRTP.h"

struct I::BeamformerAdaptive
{
	Complex2D Input;
	Boolean SpeechActivity;
};

class BeamformerAdaptive : public Base<BeamformerAdaptive, I::BeamformerAdaptive, O::Complex>
{
	friend Base<BeamformerAdaptive, I::BeamformerAdaptive, O::Complex>;
protected:
	void CalculateFilter()
	{
		D.Rx[D.CurrentBand] += 1e-16f*Eigen::MatrixXcf::Identity(C.NChannels, C.NChannels);
		D.Rn[D.CurrentBand] += 1e-17f*Eigen::MatrixXcf::Identity(C.NChannels, C.NChannels);
		D.EigenSolver.compute(D.Rx[D.CurrentBand], D.Rn[D.CurrentBand]);
		D.EigenVectors = D.EigenSolver.eigenvectors();
		const float px = std::real((D.EigenVectors.col(C.NChannels - 1).adjoint() * D.Rx[D.CurrentBand] * D.EigenVectors.col(C.NChannels - 1))(0));
		const float pn = std::real((D.EigenVectors.col(C.NChannels - 1).adjoint() * D.Rn[D.CurrentBand] * D.EigenVectors.col(C.NChannels - 1))(0));
		const float gain = std::max(1.f / (1.f + pn / std::max(px - pn, 1e-20f) * P.SpeechDistortionRatio), D.WienerGainMinimum); // Wiener gain written so it is stable when px-pn=0 and P.SpeechDistortionRatio=0 (should give gain=1)
		D.EigenVectors.col(C.NChannels - 1) *= -D.EigenVectors.block(1, 0, C.NChannels - 1, C.NChannels - 1).determinant() / D.EigenVectors.determinant();
		// put resulting beamformer for current band into Filter 
		D.Filter.row(D.CurrentBand) = D.EigenVectors.col(C.NChannels - 1).conjugate() * gain;
	}

private:
	struct Coefficients
	{
		int NChannels = 2;
		float FilterbankRate = 125.f;
		int NBands = 257;
	} C;

	struct Parameters
	{
		enum SpeechDecisions { Noise, Speech, Automatic, FreezeUpdate };
		SpeechDecisions SpeechDecision = Automatic;
		ConstrainedType<float> FilterUpdateTConstant = { 5.f, 0.001f, 60.f }; // seconds
		ConstrainedType<float> FilterUpdateRate = { 8.f, 0.001f, 100.f }; // Hz
		ConstrainedType<float> SpeechDistortionRatio = { 0.f, 0.f, 100.f };
		ConstrainedType<float> WienerGainMinimumdB = { -10.f, -60.f, 0.f };
	} P;

	struct Data
	{
		std::vector<Eigen::MatrixXcf> Rx, Rn;
		Eigen::ArrayXXcf Filter;
		Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcf> EigenSolver;
		Eigen::MatrixXcf EigenVectors, Rxn;
		int CurrentBand, FilterUpdatesPerFrame;
		float FilterbankRate, FilterUpdateLambda, WienerGainMinimum;
		bool ActivityFlag;

		void Reset()
		{
			CurrentBand = 0;
			ActivityFlag = false;
			EigenVectors.setZero();
			Filter.setZero();
			Filter.col(0) = 1;
			for (auto i = 0u; i < Rx.size(); i++)
			{
				Rx[i].setZero();
				Rx[i](0,0) = 1e-5f;
				Rn[i].setZero();
			}
			Rxn.setZero();
			EigenSolver.compute(Rx[0], Rn[0]); // run once to make sure calculation works
		};
		bool InitializeMemory(const Coefficients& c)
		{
			EigenVectors.resize(c.NChannels, c.NChannels);
			Filter.resize(c.NBands, c.NChannels);
			Rx.resize(c.NBands);
			for (auto& rx : Rx) { rx.resize(c.NChannels, c.NChannels); }
			Rn.resize(c.NBands);
			for (auto& rn : Rn) { rn.resize(c.NChannels, c.NChannels); }
			Rxn.resize(c.NChannels, c.NChannels);
			EigenSolver = Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcf>::GeneralizedSelfAdjointEigenSolver(c.NChannels);
			return true;
		}
		size_t GetAllocatedMemorySize() const
		{
			auto size = Filter.GetAllocatedMemorySize();
			for (auto& rx : Rx) { size += rx.GetAllocatedMemorySize(); }
			for (auto& rn : Rn) { size += rn.GetAllocatedMemorySize(); }
			size += EigenVectors.GetAllocatedMemorySize();
			size += Rxn.GetAllocatedMemorySize();
			return size;
		}
		void OnParameterChange(const Parameters& p, const Coefficients& c)
		{
			FilterUpdateLambda = 1.f - expf(-1.f / (c.FilterbankRate * p.FilterUpdateTConstant));
			FilterUpdatesPerFrame = static_cast<int>(c.NBands * p.FilterUpdateRate / c.FilterbankRate);
			WienerGainMinimum = powf(10, p.WienerGainMinimumdB * 0.05f);
		}
	} D;

	void ProcessOn(Input xFreq, Output yFreq)
	{
		D.ActivityFlag = xFreq.SpeechActivity;

		switch (P.SpeechDecision)
		{
		case Parameters::Speech: D.ActivityFlag = true; break;
		case Parameters::Noise: D.ActivityFlag = false; break;
		default: break;
		}

		if (P.SpeechDecision != Parameters::FreezeUpdate) { CovarianceUpdate(xFreq.Input); }

		for (auto i = 0;i < D.FilterUpdatesPerFrame;i++)
		{
			CalculateFilter();
			D.CurrentBand++;
			if (D.CurrentBand >= C.NBands) { D.CurrentBand = 0; }
		}

		yFreq = (xFreq.Input * D.Filter).rowwise().sum(); // this has been profiled to be just as fast as multiplying with ones or summing in a for-loop
	}

	void ProcessOff(Input xFreq, Output yFreq) { yFreq = xFreq.Input.col(0); }

	void CovarianceUpdate(I::Complex2D Input)
	{
		for (auto band = 0;band < C.NBands;band++)
		{
			// only calculate lower triangular covariance matrix, since that is all that is used by eigensolver
			// Full matrix is: D.Rxn = Input.matrix().row(i).transpose() * Input.matrix().row(i).conjugate();
			for (auto channel = 0;channel < C.NChannels;channel++)
			{
				D.Rxn.block(channel, channel, C.NChannels - channel, 1) = Input.block(band,channel,1,C.NChannels-channel).transpose().matrix() * std::conj(Input(band, channel));
			}
			D.Rx[band] += D.FilterUpdateLambda * (D.Rxn - D.Rx[band]);
			if (!D.ActivityFlag) { D.Rn[band] += D.FilterUpdateLambda * (D.Rxn - D.Rn[band]); }
		}
	}
};