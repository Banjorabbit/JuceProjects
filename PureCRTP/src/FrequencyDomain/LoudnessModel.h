#pragma once
#include "../BaseClasses/PureCRTP.h"
#include "CriticalBands.h"

// LoudnessModel is a simple loudness model that can be used to calculate a more perceptual representation of activity in the signal.
//
// author: Kristian Timm Andersen
class LoudnessModel : public Base<LoudnessModel>
{
	friend Base<LoudnessModel>;

public:
	CriticalBands ConvertBands;

private:

	struct Coefficients 
	{
		float SampleRate = 44.1e3f;
		float FilterbankRate = 44.1e3f / 256;
		int NBands = 257;
		int NChannels = 2;
	} C;

	struct Parameters 
	{
		float AttackTConstant = .01f;
		float ReleaseTConstant = .05f;
	} P;

	struct Data 
	{
		Eigen::ArrayXf Weighting; // frequency weighting
		int NBandsCritical; 
		float AttackLambda, ReleaseLambda;
		Eigen::ArrayXXf LoudnessPeak;
		void Reset() { LoudnessPeak.setZero(); }
		bool InitializeMemory(const Coefficients& c)
		{
			// calculate frequency weighting
			Weighting.resize(c.NBands);
			auto NFFT = (c.NBands - 1) * 2; // number of FFT points
			float FreqResolution = c.SampleRate / NFFT;
			auto n150 = static_cast<int>(std::ceilf(150.f / FreqResolution));
			Weighting.head(n150) = Eigen::ArrayXf::LinSpaced(n150, 0.f, 1.f);
			auto n750 = static_cast<int>(750.f / FreqResolution);
			Weighting.segment(n150, n750 - n150).setOnes();
			auto n2500 = static_cast<int>(2500.f / FreqResolution);
			Weighting.segment(n750, n2500 - n750) = Eigen::ArrayXf::LinSpaced(n2500 - n750, 1.f, 2.5f);
			Weighting.segment(n2500, c.NBands - n2500).setConstant(2.5f);
			return true;
		}
		size_t GetAllocatedMemorySize() const { return LoudnessPeak.GetAllocatedMemorySize() + Weighting.GetAllocatedMemorySize(); }
		void OnParameterChange(const Parameters& p, const Coefficients& c) 
		{
			AttackLambda = 1.f - expf(-1.f / (c.FilterbankRate*p.AttackTConstant));
			ReleaseLambda = 1.f - expf(-1.f / (c.FilterbankRate*p.ReleaseTConstant));
		}
	} D;

	DEFINEMEMBERALGORITHMS(1, ConvertBands);

	auto InitializeMembers()
	{
		auto c = ConvertBands.GetCoefficients();
		c.NBands = C.NBands;
		c.SampleRate = C.SampleRate;
		auto flag =  ConvertBands.Initialize(c);
		D.NBandsCritical = ConvertBands.GetNBandsCritical();
		D.LoudnessPeak.resize(D.NBandsCritical, C.NChannels);
		D.LoudnessPeak.setZero();
		return flag;
	}

	void ProcessOn(Input xPower, Output yLoudness) 
	{
		for (auto channel = 0; channel < C.NChannels; channel++)
		{
			Eigen::ArrayXf powerCritical(D.NBandsCritical);
			ConvertBands.Process(xPower.col(channel) * D.Weighting, powerCritical);
			for (auto i = 0; i < D.NBandsCritical; i++)
			{
				// take power x^0.3
				float loudness = powf(powerCritical(i), 0.3f);
				// 1st order lowpass with assymmetric time constants
				loudness -= D.LoudnessPeak(i,channel);
				D.LoudnessPeak(i, channel) += loudness > 0 ? D.AttackLambda * loudness : D.ReleaseLambda * loudness;
			}
		}
		yLoudness = D.LoudnessPeak;
	}

	void ProcessOff(Input xPower, Output yLoudness) { yLoudness.setZero(); }
};
