#pragma once
#include "../BaseClasses/PureCRTP.h"
#include "NoiseEstimation.h"

// Detect voice activation based on the short-time energy. The output is a single bool that integrates the voice activation accross channels and frequency bands. 
// If detection is needed for each frequency band, then use getActivityArray()
//
// author: Kristian Timm Andersen
class DetectVoiceActivation : public Base<DetectVoiceActivation, I::Complex2D, O::Boolean>
{
	friend Base<DetectVoiceActivation, I::Complex2D, O::Boolean>;

public:
	VectorAlgo<NoiseEstimationSPP> NoiseActivities;

	Eigen::Array<bool, Eigen::Dynamic, 1> getActivityArray() const { return D.activityArray; }

private:
	struct Coefficients 
	{
		float FilterbankRate = 125.f;
		int NBands = 257;
		int NChannels = 2;
	} C;

	struct Parameters
	{
		ConstrainedType<float> SnrThreshold = { 1.5f, 0.f, 100.f };
		ConstrainedType<float> ActivityThreshold = { 0.5f, 0.f, 100.f };
		ConstrainedType<float> ActivityMeanThreshold = { 0.15f, 0.f, 100.f };
		ConstrainedType<float> VADThreshold = { 1.f, 0.f, 100.f };
	} P;
	
	struct Data
	{
		Eigen::Array<bool, Eigen::Dynamic, 1> activityArray;
		void Reset() 
		{ 
			activityArray.setConstant(false);
		}
		bool InitializeMemory(const Coefficients& c) 
		{ 
			activityArray.resize(c.NBands);
			return true; 
		}
		size_t GetAllocatedMemorySize() const 
		{ 
			size_t size = 0;
			size += activityArray.GetAllocatedMemorySize();
			return size;
		}
		void OnParameterChange(const Parameters& p, const Coefficients& c)	{}
	} D;

	DEFINEMEMBERALGORITHMS(1, NoiseActivities);

	auto InitializeMembers()
	{ 
		NoiseActivities.resize(C.NChannels);
		auto cNA = NoiseActivities[0].GetCoefficients();
		cNA.FilterbankRate = C.FilterbankRate;
		cNA.NBands = C.NBands;
		cNA.NChannels = 1;
		return NoiseActivities.Initialize(cNA);
	}

	void ProcessOn(Input xFreq, Output activityFlag)
	{
		D.activityArray.setConstant(false);
		activityFlag = false;
		// fuse channels into one activity detector
		for (auto channel = 0; channel < xFreq.cols(); channel++)
		{
			Eigen::ArrayXf activityChannel(C.NBands);
			Eigen::ArrayXf powerNoisy = xFreq.col(channel).abs2();
			Eigen::ArrayXf powerNoise(C.NBands);
			NoiseActivities[channel].Process(powerNoisy, { powerNoise, activityChannel }); // Run noise estimator and activity detector
			activityFlag |= powerNoisy.mean() > (powerNoise.mean() * P.SnrThreshold * P.VADThreshold);
			D.activityArray = D.activityArray.operator||(activityChannel > (P.ActivityThreshold * P.VADThreshold));
		}
		activityFlag |= static_cast<float>(D.activityArray.count()) > (P.ActivityMeanThreshold * P.VADThreshold * C.NBands);
	}

	void ProcessOff(Input xFreq, Output activityFlag) { activityFlag = false; }
};