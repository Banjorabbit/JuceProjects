#pragma once
#include "../BaseClasses/PureCRTP.h"
#include "NoiseEstimation.h"
#include "GainCalculation.h"

class NoiseSuppression : public BaseFrequencyDomain<NoiseSuppression>
{
	friend BaseFrequencyDomain<NoiseSuppression>;

public:
	NoiseEstimationSPP noiseEstimation;
	GainCalculationSimple gainCalculation;

private:
	struct Coefficients	
	{
		float FilterbankRate = 125;
		int NBands = 257;
		int NChannels = 1;
	} C;

	struct Parameters {} P;

	struct Data
	{
		void Reset() { }
		bool InitializeMemory(const Coefficients& c) { return true; }
		size_t GetAllocatedMemorySize() const { return 0; }
		void OnParameterChange(const Parameters& p, const Coefficients& c) {}
	} D;

	DEFINEMEMBERALGORITHMS(2, noiseEstimation, gainCalculation);

	auto InitializeMembers() 
	{
		auto cNE = noiseEstimation.GetCoefficients();
		cNE.FilterbankRate = C.FilterbankRate;
		cNE.NBands = C.NBands;
		cNE.NChannels = C.NChannels;
		auto flag = noiseEstimation.Initialize(cNE);

		auto cGC = gainCalculation.GetCoefficients();
		cGC.FilterbankRate = C.FilterbankRate;
		cGC.NBands = C.NBands;
		cGC.NChannels = C.NChannels;
		flag &= gainCalculation.Initialize(cGC);
		return flag;
	}

	void ProcessOn(Input xFreq, Output yFreq)
	{
		Eigen::ArrayXXf power = xFreq.abs2();
		Eigen::ArrayXXf powerNoise(C.NBands, C.NChannels);
		Eigen::ArrayXXf gain(C.NBands, C.NChannels);
		noiseEstimation.Process(power, { powerNoise, gain }); // second output is not used, so write to gain that will be overwritten later
		power /= (powerNoise + 1e-20f);
		gainCalculation.Process(power, gain);
		yFreq = xFreq * gain;
	}

	void ProcessOff(Input xFreq, Output yFreq)
	{
		yFreq = xFreq;
	}
};