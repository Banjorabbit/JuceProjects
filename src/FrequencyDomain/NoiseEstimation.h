#pragma once
#include "../BaseClasses/PureCRTP.h"

struct O::NoiseEstimationSPP
{
	Real2D PowerNoise;
	Real2D Activity;
};

class NoiseEstimationSPP : public Base<NoiseEstimationSPP, I::Real2D, O::NoiseEstimationSPP>
{
	friend Base<NoiseEstimationSPP, I::Real2D, O::NoiseEstimationSPP>;

	struct Coefficients
	{
		float FilterbankRate = 125.f;
		int NBands = 257;
		int NChannels = 2;
		float SnrAprioridB = 15;
		float SpeechPresenceApriori = 0.5f;
	} C;

	struct Parameters
	{
		ConstrainedType<float> ActivityMeanTConstant = { .152f, 0.001f, 60.f }; // seconds
		ConstrainedType<float> NoiseUpdateTConstant = { .072f, 0.001f, 60.f }; // seconds
	} P;

	struct Data
	{
		Eigen::ArrayXXf ActivityMean;
		float ActivityMeanLambda, NoiseUpdateLambda, Offset, Scaling, ScalingLin;
		void Reset() { ActivityMean.setConstant(0.5f); };
		bool InitializeMemory(const Coefficients& c)
		{
			ActivityMean.resize(c.NBands, c.NChannels);
			return true;
		}
		size_t GetAllocatedMemorySize() const { return ActivityMean.GetAllocatedMemorySize(); }
		void OnParameterChange(const Parameters& p, const Coefficients& c)
		{
			ActivityMeanLambda = 1.f - expf(-1.f / (c.FilterbankRate*p.ActivityMeanTConstant));
			NoiseUpdateLambda = 1.f - expf(-1.f / (c.FilterbankRate*p.NoiseUpdateTConstant));
			const float snrApriori = powf(10, c.SnrAprioridB * .1f);
			Offset = log(1.f / (1.f + snrApriori));
			Scaling = snrApriori / (1.f + snrApriori);
			ScalingLin = c.SpeechPresenceApriori / (1.f - c.SpeechPresenceApriori);
		}
	} D;

	void ProcessOn(Input powerNoisy, Output y)
	{
		// Calculate activity detection
		y.Activity = D.ScalingLin * (D.Scaling * powerNoisy / y.PowerNoise + D.Offset).cwiseMin(25.f).exp();
		y.Activity /= (1.f + y.Activity);

		// Update long term activity
		D.ActivityMean += D.ActivityMeanLambda * (y.Activity - D.ActivityMean);
		y.Activity = (D.ActivityMean > 0.99f).select(y.Activity.cwiseMin(.99f), y.Activity);

		// Update Noise
		y.PowerNoise += D.NoiseUpdateLambda * (1.f - y.Activity) * (powerNoisy - y.PowerNoise);
	}

	void ProcessOff(Input powerNoisy, Output y) { y.PowerNoise = powerNoisy; y.Activity.setZero(); }

};

