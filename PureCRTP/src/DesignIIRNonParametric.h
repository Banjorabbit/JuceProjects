#pragma once
#include "BaseClasses/PureCRTP.h"
#include "CubicSpline.h"
#include "DesignIIRMinPhase.h"

struct I::DesignIIRNonParametric
{
	Real2D Frequencies;
	Real2D GaindB;

};

struct O::DesignIIRNonParametric
{
	Real2D SOS;
	Float Gain;
};

// Calculate minimum phase IIR filter from a number of frequency points and corresponding gain in dB.
//
// author: Kristian Timm Andersen
class DesignIIRNonParametric : public Base<DesignIIRNonParametric, I::DesignIIRNonParametric, O::DesignIIRNonParametric>
{
	friend Base<DesignIIRNonParametric, I::DesignIIRNonParametric, O::DesignIIRNonParametric>;

public:
	CubicSpline SplineCalculation;
	DesignIIRMinPhase FilterDesigner;

private:
	struct Coefficients
	{
		int N = 8;
		int NBands = 65;
		float SampleRate = 16e3f;
	} C;

	struct Parameters { } P;

	struct Data
	{
		void Reset() {}
		bool InitializeMemory(const Coefficients& c) { return true; }
		size_t GetAllocatedMemorySize() const { return 0; }
		void OnParameterChange(const Parameters& p, const Coefficients& c) { }
	} D;

	DEFINEMEMBERALGORITHMS(2, SplineCalculation, FilterDesigner);

	auto InitializeMembers()
	{
		// spline
		auto sSpline = SplineCalculation.GetSetup();
		sSpline.Parameters.sfactor = 0.f; // set spline interpolation to Catmull-Rom spline
		auto flag = SplineCalculation.Initialize(sSpline);
		// min phase filter designer
		auto setupFD = FilterDesigner.GetSetup();
		setupFD.Coefficients.N = C.N;
		setupFD.Coefficients.WeightType = setupFD.Coefficients.MelScale;
		setupFD.Coefficients.NBands = C.NBands;
		setupFD.Coefficients.SampleRate = C.SampleRate;
		flag &= FilterDesigner.Initialize(setupFD);
		return flag;
	}

	void ProcessOn(Input xFreq, Output yTime)
	{
		Eigen::ArrayXf FreqsFFT = Eigen::ArrayXf::LinSpaced(C.NBands, 0.f, C.SampleRate / 2); // desired frequency points for FFT
		// for each filter
		for (auto channel = 0; channel < xFreq.Frequencies.cols(); channel++)
		{
			// add 2 frequency/gain points to the left and right of given points to ensure correct gradient 
			Eigen::ArrayXf Frequencies(xFreq.Frequencies.rows() + 4);
			Frequencies.head(2) = -xFreq.Frequencies.col(channel).head(2).reverse();
			Frequencies.segment(2, xFreq.Frequencies.rows()) = xFreq.Frequencies.col(channel);
			Frequencies.tail(2) = C.SampleRate - xFreq.Frequencies.col(channel).tail(2).reverse();
			Eigen::ArrayXf GaindB(xFreq.GaindB.rows() + 4);
			GaindB.head(2) = xFreq.GaindB.col(channel).head(2).reverse();
			GaindB.segment(2, xFreq.GaindB.rows()) = xFreq.GaindB.col(channel);
			GaindB.tail(2) = xFreq.GaindB.col(channel).tail(2).reverse();
			// calculate interpolation
			Eigen::ArrayXf GainFFT(C.NBands);
			SplineCalculation.Process({ Frequencies, GaindB, FreqsFFT }, GainFFT);
			// convert to linear scale
			GainFFT *= 0.115129254649702f; // 10^(GainFFT/20) = e^(log(10^(GainFFT/20))) = e^(GainFFT/20*log(10)) = e^(GainFFT*0.115129254649702)
			GainFFT = GainFFT.exp();
			// calculate minimum phase spectrum
			FilterDesigner.Process(GainFFT, { yTime.SOS, yTime.Gain });
		}
	}

	void ProcessOff(Input xFreq, Output yTime) { yTime.SOS.setZero(); yTime.SOS.col(3).setOnes(); yTime.Gain = 1.f; }
};