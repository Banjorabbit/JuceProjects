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

private:
	struct Coefficients
	{
		float FilterbankRate = 125.f;
		float FilterLengthMilliseconds = 250.f;
		int NBands = 257;
		int NChannels = 2;
		int nBuffersInAnalysisFrame = 4;
		int nBuffersInSynthesisFrame = 4;
	} C;

	struct Parameters
	{
		ConstrainedType<float> NearendLimitdB = { -40.f, -100.f, 0.f }; // dB
	} P;

	struct Data
	{
		std::vector<Eigen::ArrayXXcf> Filters;
		int BufferLength;
		float Lambda, NearendLimit;
		std::vector<Eigen::ArrayXf> CoefficientVariance, Momentums;
		Eigen::ArrayXXcf BuffersLoopback;
		Eigen::ArrayXf LoopbackVariance;
		Eigen::ArrayXXf NearendVariance;
		int CircCounter;

		void Reset()
		{
			BuffersLoopback.setZero();
			for (auto &filter : Filters) { filter.setZero(); }
			LoopbackVariance.setConstant(100.f);
			for (auto &coefficientVariance : CoefficientVariance) { coefficientVariance.setZero(); }
			for (auto &momentum : Momentums) { momentum.setOnes(); }
			NearendVariance.setConstant(10.f);
			CircCounter = 0;
		}
		bool InitializeMemory(const Coefficients& c)
		{
			BufferLength = std::max(static_cast<int>(c.FilterLengthMilliseconds * c.FilterbankRate * 1e-3f) - (c.nBuffersInAnalysisFrame + c.nBuffersInSynthesisFrame - 1), 0) + 1;
			Lambda = 1.f - expf(-1.f / (c.nBuffersInAnalysisFrame + BufferLength - 1));
			BuffersLoopback.resize(BufferLength, c.NBands);
			Filters.resize(c.NChannels);
			for (auto &filter : Filters) { filter.resize(BufferLength, c.NBands); }
			LoopbackVariance.resize(c.NBands);
			CoefficientVariance.resize(c.NChannels);
			for (auto& coefficientVariance : CoefficientVariance) { coefficientVariance.resize(c.NBands); }
			Momentums.resize(c.NChannels);
			for (auto& momentum : Momentums) { momentum.resize(c.NBands); }
			NearendVariance.resize(c.NBands, c.NChannels);
			return true;
		}
		size_t GetAllocatedMemorySize() const
		{
			size_t size = 0;
			for (auto& filter : Filters) { size += filter.GetAllocatedMemorySize(); }
			size += BuffersLoopback.GetAllocatedMemorySize();
			size += LoopbackVariance.GetAllocatedMemorySize();
			for (auto& coef : CoefficientVariance) { size += coef.GetAllocatedMemorySize(); }
			for (auto& momentum : Momentums) { size += momentum.GetAllocatedMemorySize(); }
			size += NearendVariance.GetAllocatedMemorySize();
			return size;
		}
		void OnParameterChange(const Parameters& p, const Coefficients& c)
		{
			NearendLimit = powf(10.f, p.NearendLimitdB * 0.1f);
		}
	} D;

	void ProcessOn(Input xFreq, Output yFreq)
	{
		D.CircCounter = D.CircCounter <= 0 ? D.BufferLength - 1 : D.CircCounter - 1;
		// Put loopback signals into buffer
		D.BuffersLoopback.row(D.CircCounter) = xFreq.Loopback.transpose();
		
		// update loopback variance
		D.LoopbackVariance += D.Lambda * (xFreq.Loopback.abs2() - D.LoopbackVariance); // variance of loopback signal i

		//  filter and subtract from input
		const int conv_length1 = D.BufferLength - D.CircCounter;
		for (auto channel = 0; channel < C.NChannels; channel++)
		{
			// estimated transferfunction x loopback
			for (auto ibin = 0; ibin < C.NBands; ibin++)
			{
				std::complex<float> micEst = 0.f;
				for (auto i = 0; i < conv_length1; i++)
				{
					micEst += D.BuffersLoopback(D.CircCounter + i, ibin) * D.Filters[channel](i, ibin);
				}
				for (auto i = 0; i < D.CircCounter; i++)
				{
					micEst += D.BuffersLoopback(i, ibin) * D.Filters[channel](conv_length1 + i, ibin);
				}
				yFreq(ibin, channel) = xFreq.Input(ibin, channel) - micEst;

				const float newpow = (micEst.real()*micEst.real() + micEst.imag()*micEst.imag()) * D.NearendLimit;
				const float pow_error = (yFreq(ibin, channel).real()*yFreq(ibin, channel).real() + yFreq(ibin, channel).imag()*yFreq(ibin, channel).imag());
				D.NearendVariance(ibin, channel) += D.Lambda * (std::max(pow_error, newpow) - D.NearendVariance(ibin, channel));

				const float p = D.Momentums[channel](ibin) + D.BufferLength * D.CoefficientVariance[channel](ibin);
				float q = p / (D.BufferLength * D.NearendVariance(ibin, channel) + (D.BufferLength + 2) * p * D.LoopbackVariance(ibin) + 1e-30f);
				q = std::min(q, 1.f / (D.LoopbackVariance(ibin) * D.BufferLength + 1e-30f));
				D.Momentums[channel](ibin) = std::max((1.f - q * D.LoopbackVariance(ibin)) * p, 0.01f);
				const std::complex<float> W = q * yFreq(ibin, channel);

				const std::complex<float> temp1 = W * std::conj(D.BuffersLoopback(D.CircCounter, ibin));
				D.Filters[channel](0, ibin) += temp1;
				for (auto i = 1; i < conv_length1; i++)
				{
					D.Filters[channel](i, ibin) += W * std::conj(D.BuffersLoopback(D.CircCounter + i, ibin));
				}
				for (auto i = 0; i < D.CircCounter; i++)
				{
					D.Filters[channel](conv_length1 + i, ibin) += W * std::conj(D.BuffersLoopback(i, ibin));
				}
				// update coefficient variance
				D.CoefficientVariance[channel](ibin) += D.Lambda * (temp1.real()*temp1.real() + temp1.imag()*temp1.imag() - D.CoefficientVariance[channel](ibin));
			}
		}
	}

	void ProcessOff(Input xFreq, Output yFreq) { yFreq = xFreq.Input; }
};