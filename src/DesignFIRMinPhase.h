#pragma once
#include "BaseClasses/PureCRTP.h"
#include "FrequencyDomain/MinPhaseSpectrum.h"
#include "FFT.h"

// Calculate minimum phase FIR filter a real gain vector.
//
// author: Kristian Timm Andersen
class DesignFIRMinPhase : public Base<DesignFIRMinPhase>
{
	friend Base<DesignFIRMinPhase>;

public:
	MinPhaseSpectrum MinPhaseCalculation;
	FFTReal FFT;

private:
	struct Coefficients
	{
		int FilterSize = 128;
	} C;

	struct Parameters {} P;

	struct Data
	{
		int NBands;
		void Reset() {}
		bool InitializeMemory(const Coefficients& c)
		{
			NBands = c.FilterSize / 2 + 1;
			return true;
		}
		size_t GetAllocatedMemorySize() const { return 0; }
		void OnParameterChange(const Parameters& p, const Coefficients& c) { }
	} D;

	DEFINEMEMBERALGORITHMS(2, MinPhaseCalculation, FFT);

	auto InitializeMembers()
	{
		// min phase calculation
		auto setupMPS = MinPhaseCalculation.GetSetup();
		setupMPS.Coefficients.NBands = D.NBands;
		auto flag = MinPhaseCalculation.Initialize(setupMPS);
		// FFT
		auto setupFFT = FFT.GetSetup();
		setupFFT.Coefficients.FFTSize = C.FilterSize;
		flag &= FFT.Initialize(setupFFT);
		return flag;
	}

	void ProcessOn(Input xMag, Output yTime)
	{
		// for each filter
		for (auto channel = 0; channel < xMag.cols(); channel++)
		{
			// calculate minimum phase spectrum
			Eigen::ArrayXcf MinPhaseSpectrum(D.NBands);
			MinPhaseCalculation.Process(xMag.col(channel), MinPhaseSpectrum);
			// convert to time domain filter
			FFT.Inverse(MinPhaseSpectrum, yTime.col(channel));
		}
	}

	void ProcessOff(Input xMag, Output yTime) { yTime.setZero(); }
};