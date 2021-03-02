#pragma once
#include "../BaseClasses/PureCRTP.h"
#include "CriticalBands.h"

// Detect transients in critical bands. The output mapped to critical bands, 
// but can be mapped back to the original bands using 
// DetectTransient.ConvertBands.Inverse()
//
// author: Kristian Timm Andersen
class DetectTransient : public Base<DetectTransient, I::Real2D, O::Bool2D>
{
	friend Base<DetectTransient, I::Real2D, O::Bool2D>;

public:
	CriticalBands ConvertBands;

private:
	struct Coefficients 
	{
		int NBands = 257;
		float SampleRate = 44.1e3f;
		float FilterbankRate = 44.1e3f / 256;
		int NChannels = 2;
	} C;

	struct Parameters 
	{
		float AttackShortTConstant = .001f;
		float ReleaseShortTConstant = .01f;
		float LongTConstant = .025f;
		ConstrainedType<float> TransientThreshold = { 0.75f, 0.f, 1.f };
	} P;

	struct Data 
	{
		int NBandsCritical;
		Eigen::ArrayXXf PerceptualShort;
		Eigen::ArrayXXf PerceptualLong;
		float AttackShortLambda, ReleaseShortLambda, LongLambda;
		void Reset() 
		{
			PerceptualShort.setZero();
			PerceptualLong.setZero();
		}
		bool InitializeMemory(const Coefficients& c)
		{
			return true;
		}
		size_t GetAllocatedMemorySize() const 
		{ 
			size_t size = PerceptualShort.GetAllocatedMemorySize();
			size += PerceptualLong.GetAllocatedMemorySize();
			return size; 
		}
		void OnParameterChange(const Parameters& p, const Coefficients& c) 
		{
			AttackShortLambda = 1.f - expf(-1.f / std::max(1e-10f, c.FilterbankRate * p.AttackShortTConstant));
			ReleaseShortLambda = 1.f - expf(-1.f / std::max(1e-10f, c.FilterbankRate * p.ReleaseShortTConstant));
			LongLambda = 1.f - expf(-1.f / std::max(1e-10f, c.FilterbankRate * p.LongTConstant));
		}
	} D;

	DEFINEMEMBERALGORITHMS(1, ConvertBands);

	auto InitializeMembers()
	{
		auto c = ConvertBands.GetCoefficients();
		c.NBands = C.NBands;
		c.SampleRate = C.SampleRate;
		bool flag = ConvertBands.Initialize(c);
		D.NBandsCritical = ConvertBands.GetNBandsCritical();
		D.PerceptualShort.resize(C.NBands, C.NChannels);
		D.PerceptualShort.setZero();
		D.PerceptualLong.resize(C.NBands, C.NChannels);
		D.PerceptualLong.setZero();
		return flag;
	}

	void ProcessOn(Input xPower, Output activityFlag) 
	{
		for (auto channel = 0; channel < xPower.cols(); channel++)
		{
			Eigen::ArrayXf xPerceptual(D.NBandsCritical);
			ConvertBands.Process(xPower.col(channel), xPerceptual);
			for (auto i = 0; i < ConvertBands.GetNBandsCritical(); i++)
			{
				auto xShortDiff = xPerceptual(i) - D.PerceptualShort(i);
				D.PerceptualShort(i) += xShortDiff > 0.f ? D.AttackShortLambda * xShortDiff : D.ReleaseShortLambda * xShortDiff;
				D.PerceptualLong(i) += D.LongLambda * (xPerceptual(i) - D.PerceptualLong(i));
				auto activityProbability = D.PerceptualShort(i) / (D.PerceptualShort(i) + D.PerceptualLong(i));
				activityFlag(i, channel) = activityProbability > P.TransientThreshold;
			}
		}
	}

	void ProcessOff(Input xPower, Output activityFlag) { activityFlag.setConstant(0); }
};