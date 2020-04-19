#pragma once
#include "../BaseClasses/PureCRTP.h"
#include "NoiseEstimation.h"
#include "GainCalculation.h"

class NoiseSuppression : public BaseFrequencyDomain<NoiseSuppression>
{
	friend BaseFrequencyDomain<NoiseSuppression>;

public:
	VectorAlgo<NoiseEstimationSPP> noiseEstimation;
	VectorAlgo<GainCalculationSimple> gainCalculation;

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
		Eigen::ArrayXXf PowerNoise;
		Eigen::ArrayXXf Gain;
		void Reset() { Gain.setZero(); PowerNoise.setZero(); }
		bool InitializeMemory(const Coefficients& c) 
		{
			Gain.resize(c.NBands, c.NChannels);
			PowerNoise.resize(c.NBands, c.NChannels);
			return true;
		}
		size_t GetAllocatedMemorySize() const { return PowerNoise.GetAllocatedMemorySize() + Gain.GetAllocatedMemorySize(); }
		void OnParameterChange(const Parameters& p, const Coefficients& c) {}
	} D;

	DEFINEMEMBERALGORITHMS(2, noiseEstimation, gainCalculation);

	auto InitializeMembers() 
	{
		noiseEstimation.resize(C.NChannels);
		gainCalculation.resize(C.NChannels);
		auto cNE = noiseEstimation[0].GetCoefficients();
		auto cGC = gainCalculation[0].GetCoefficients();
		cNE.FilterbankRate = C.FilterbankRate;
		cNE.NBands = C.NBands;
		cNE.NChannels = 1;
		cGC.FilterbankRate = C.FilterbankRate;
		auto flag = noiseEstimation.Initialize(cNE);
		flag &= gainCalculation.Initialize(cGC);
		return flag;
	}

	void ProcessOn(Input xFreq, Output yFreq)
	{
		for (auto channel = 0; channel < xFreq.cols(); channel++)
		{
			Eigen::ArrayXf power = xFreq.col(channel).abs2();
			Eigen::ArrayXf activity(C.NBands);
			noiseEstimation[channel].Process(power, { D.PowerNoise.col(channel), activity });
			power /= D.PowerNoise.col(channel);
			gainCalculation[channel].Process(power, D.Gain.col(channel));
		}
		yFreq = xFreq * D.Gain;
	}

	void ProcessOff(Input xFreq, Output yFreq)
	{
		yFreq = xFreq;
	}
};