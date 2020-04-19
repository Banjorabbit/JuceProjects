#pragma once
#include "../BaseClasses/PureCRTP.h"

class GainCalculationSimple : public Base<GainCalculationSimple>
{
	friend Base<GainCalculationSimple>;

	struct Coefficients
	{
		float FilterbankRate = 125;
	} C;

	struct Parameters
	{
		ConstrainedType<float> AprioriSmoother = { .98f, 0.f, 1.f };
		ConstrainedType<float> SpeechDistortionRatio = { 1.f, 0.f, 100.f };
		ConstrainedType<float> MinimumdB = { -15.f, -100.f, 0.f }; // dB
		ConstrainedType<float> Exponential = { 1.4f, 0.f, 10.f };
		ConstrainedType<float> AttackTConstant = { .001f, 1e-6f, 60.f }; // seconds
		ConstrainedType<float> ReleaseTConstant = { 0.02f, 1e-6f, 60.f }; // seconds
		enum Methods { SimpleApriori, SimpleTwoStep, Smoothing};
		Methods Method = SimpleApriori;
	} P;

	struct Data
	{
		float Minimum, AttackLambda, ReleaseLambda;
		void Reset() {}
		bool InitializeMemory(const Coefficients& c) { return true; }
		size_t GetAllocatedMemorySize() const { return 0; }
		void OnParameterChange(const Parameters& p, const Coefficients& c) 
		{
			Minimum = powf(10.f, p.MinimumdB * .05f);
			AttackLambda = 1.f - expf(-1.f / (c.FilterbankRate*p.AttackTConstant));
			ReleaseLambda = 1.f - expf(-1.f / (c.FilterbankRate*p.ReleaseTConstant));
		}
	} D;

	void ProcessOn(Input snrAposteriori, Output gain)
	{
		Eigen::ArrayXXf power = snrAposteriori.min(100.f).max(1.01f);
		switch (P.Method)
		{
		case Parameters::Methods::SimpleApriori :
			power *= (P.AprioriSmoother * (gain*gain - 1.f) + 1.f);
			power += P.AprioriSmoother - 1.f; // snr
			break;
		case Parameters::Methods::SimpleTwoStep :
			power *= (P.AprioriSmoother * (gain*gain - 1.f) + 1.f);
			power += P.AprioriSmoother - 1.f; // snr
			power = (power / (power + P.SpeechDistortionRatio)).pow(P.Exponential).max(D.Minimum); // gain
			power *= snrAposteriori * power; // two step snr
			break;
		case Parameters::Methods::Smoothing :
			power -= 1.f; // snr
			break;
		default :
			power -= 1.f; // default is the same as smoothing
		}
		power = (power / (power + P.SpeechDistortionRatio)).pow(P.Exponential).max(D.Minimum); // gain
		power -= gain; // gain difference
		gain += (power > 0).select(D.AttackLambda*power, D.ReleaseLambda*power); // smooth gain in time
	}

	void ProcessOff(Input snrAposteriori, Output gain) { gain.setOnes(); }
};
