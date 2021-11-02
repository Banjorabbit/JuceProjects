#pragma once
#include "../BaseClasses/PureCRTP.h"

struct I::EchoCancellerNLMSMinPow
{
	Complex2D Input;
	Complex Loopback;
};

class EchoCancellerNLMSMinPow : public BaseFrequencyDomain<EchoCancellerNLMSMinPow, I::EchoCancellerNLMSMinPow>
{
	friend BaseFrequencyDomain<EchoCancellerNLMSMinPow, I::EchoCancellerNLMSMinPow>;


private:
	struct Coefficients
	{
		float FilterbankRate = 125.f;
		int NBands = 257;
		int NChannels = 4;
		int FilterLength[5] = { 2, 4, 8, 16, 32 };
		ConstrainedType<int> NFilters = { 5, 1, sizeof(FilterLength) / sizeof(FilterLength[0]) }; // number of filters can not be higher than size of filterlengths
	} C;

	struct Parameters
	{
		ConstrainedType<float> NearendLimitdB = { -40.f, -100.f, 0.f }; // dB
	} P;

	struct Data
	{
		Eigen::ArrayXXcf BuffersLoopback;
		Eigen::ArrayXf Lambda;
		Eigen::ArrayXXf LoopbackVariance;
		std::vector<std::vector<Eigen::ArrayXXcf>> Filters;
		float NearendLimit;
		std::vector<Eigen::ArrayXXf> NearendVariance;
		std::vector<Eigen::ArrayXXf> Momentums;
		std::vector<Eigen::ArrayXXf> CoefficientVariance;
		int CircCounter;

		void Reset()
		{
			BuffersLoopback.setZero();
			LoopbackVariance.setConstant(1.f);
			for (auto &filter : Filters) { for (auto &f : filter) { f.setZero(); } }
			for (auto &neVar : NearendVariance) { neVar.setConstant(10.f); }
			for (auto &mom : Momentums) { mom.setOnes(); }
			for (auto &coVar : CoefficientVariance) { coVar.setZero(); }

			CircCounter = 0;
		}
		bool InitializeMemory(const Coefficients& c)
		{
			int filterLengthMax = 0;
			for (auto i = 0; i < c.NFilters; i++)
			{
				filterLengthMax = std::max(filterLengthMax, c.FilterLength[i]);
			}

			BuffersLoopback.resize(filterLengthMax, c.NBands);

			Lambda.resize(c.NFilters);
			for (auto i = 0; i < c.NFilters; i++)
			{
				Lambda(i) = 1.f - expf(-1.f / (c.FilterbankRate*0.1f + c.FilterLength[i] - 1));
			}

			LoopbackVariance.resize(c.NFilters, c.NBands);

			Filters.resize(c.NChannels);
			for (auto &filter : Filters)
			{
				filter.resize(c.NFilters);
				for (auto ifilt = 0; ifilt < c.NFilters; ifilt++)
				{
					filter[ifilt].resize(c.FilterLength[ifilt], c.NBands);
				}
			}

			NearendVariance.resize(c.NChannels);
			for (auto &neVar : NearendVariance) { neVar.resize(c.NFilters, c.NBands); }

			Momentums.resize(c.NChannels);
			for (auto &mom : Momentums) { mom.resize(c.NFilters, c.NBands); }

			CoefficientVariance.resize(c.NChannels);
			for (auto &coVar : CoefficientVariance) { coVar.resize(c.NFilters, c.NBands); }

			return true;
		}
		size_t GetAllocatedMemorySize() const
		{
			size_t size = 0;

			size += BuffersLoopback.GetAllocatedMemorySize();
			size += Lambda.GetAllocatedMemorySize();
			size += LoopbackVariance.GetAllocatedMemorySize();

			for (auto& filter : Filters) { for (auto &f : filter) { size += f.GetAllocatedMemorySize(); } }
			for (auto& neVar : NearendVariance) { size += neVar.GetAllocatedMemorySize(); }
			for (auto& mom : Momentums) { size += mom.GetAllocatedMemorySize(); }
			for (auto& coVar : CoefficientVariance) { size += coVar.GetAllocatedMemorySize(); }

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
