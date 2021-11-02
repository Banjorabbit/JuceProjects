#pragma once
#include "../BaseClasses/PureCRTP.h"

struct I::EchoCancellerNLMS
{
	Complex2D Input;
	Complex Loopback;
};

class EchoCancellerNLMS : public BaseFrequencyDomain<EchoCancellerNLMS, I::EchoCancellerNLMS>
{
	friend BaseFrequencyDomain<EchoCancellerNLMS, I::EchoCancellerNLMS>;

public:

	std::vector<Eigen::ArrayXXcf> GetFilters() const { return D.Filters; }
	void SetFilters(const std::vector<Eigen::ArrayXXcf>& filters) { D.Filters = filters; }

private:
	struct Coefficients
	{
		float FilterbankRate = 125.f;
		int FilterLength = 100;
		int NBands = 257;
		int NChannels = 4;
	} C;

	struct Parameters
	{
		ConstrainedType<float> NearendLimitdB = { -40.f, -100.f, 0.f }; // dB
	} P;

	struct Data
	{
		std::vector<Eigen::ArrayXXcf> Filters;
		float Lambda, NearendLimit;
		Eigen::ArrayXXcf BuffersLoopback;
		Eigen::ArrayXf LoopbackVariance;
		Eigen::ArrayXXf NearendVariance, Momentums, CoefficientVariance;
		int CircCounter;

		void Reset()
		{
			BuffersLoopback.setZero();
			for (auto &filter : Filters) { filter.setZero(); }
			LoopbackVariance.setConstant(100.f);
			CoefficientVariance.setZero();
			Momentums.setOnes();
			NearendVariance.setConstant(10.f);
			CircCounter = 0;
		}
		bool InitializeMemory(const Coefficients& c)
		{
			Lambda = 1.f - expf(-1.f / (c.FilterbankRate*0.1f + c.FilterLength - 1));
			BuffersLoopback.resize(c.FilterLength, c.NBands);
			Filters.resize(c.NChannels);
			for (auto &filter : Filters) { filter.resize(c.FilterLength, c.NBands); }
			LoopbackVariance.resize(c.NBands);
			CoefficientVariance.resize(c.NBands, c.NChannels);
			Momentums.resize(c.NBands, c.NChannels);
			NearendVariance.resize(c.NBands, c.NChannels);
			return true;
		}
		size_t GetAllocatedMemorySize() const
		{
			size_t size = 0;
			for (auto& filter : Filters) { size += filter.GetAllocatedMemorySize(); }
			size += BuffersLoopback.GetAllocatedMemorySize();
			size += LoopbackVariance.GetAllocatedMemorySize();
			size += CoefficientVariance.GetAllocatedMemorySize();
			size += Momentums.GetAllocatedMemorySize(); 
			size += NearendVariance.GetAllocatedMemorySize();
			return size;
		}
		void OnParameterChange(const Parameters& p, const Coefficients& c)
		{
			NearendLimit = powf(10.f, p.NearendLimitdB * 0.1f);
		}
	} D;

	void ProcessOn(Input xFreq, Output yFreq);
	void ProcessOff(Input xFreq, Output yFreq);
};
