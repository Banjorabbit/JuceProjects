#pragma once
#include "../BaseClasses/PureCRTP.h"
#include "NoiseEstimation.h"

class VoiceActivationDetection : public Base<VoiceActivationDetection, I::Complex2D, O::Boolean>
{
	friend Base<VoiceActivationDetection, I::Complex2D, O::Boolean>;

public:
	VectorAlgo<NoiseEstimationSPP> NoiseActivities;

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
		Eigen::ArrayXXf PowerNoises;
		void Reset() { PowerNoises.setConstant(1e-5f); }
		bool InitializeMemory(const Coefficients& c) { PowerNoises.resize(c.NBands, c.NChannels); return true; }
		size_t GetAllocatedMemorySize() const { return PowerNoises.GetAllocatedMemorySize(); }
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
		Eigen::Array<bool, Eigen::Dynamic, 1> activityArray = Eigen::Array<bool, Eigen::Dynamic, 1>::Constant(C.NBands, false);
		activityFlag = false;
		// fuse channels into one activity detector
		for (auto channel = 0; channel < xFreq.cols(); channel++)
		{
			Eigen::ArrayXf activityChannel(C.NBands);
			Eigen::ArrayXf powerNoisy = xFreq.col(channel).abs2();
			NoiseActivities[channel].Process(powerNoisy, { D.PowerNoises.col(channel), activityChannel }); // Run noise estimator and activity detector
			activityFlag |= powerNoisy.mean() > (D.PowerNoises.col(channel).mean() * P.SnrThreshold * P.VADThreshold);
			activityArray = activityArray.operator||(activityChannel > (P.ActivityThreshold * P.VADThreshold));
		}
		activityFlag |= static_cast<float>(activityArray.count()) > (P.ActivityMeanThreshold * P.VADThreshold * C.NBands);
	}

	void ProcessOff(Input xFreq, Output activityFlag) { activityFlag = false; }
};