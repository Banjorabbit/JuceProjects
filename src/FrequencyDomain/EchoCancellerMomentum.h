#pragma once
#include "../BaseClasses/PureCRTP.h"
#include "EchoSuppressionCovariance.h"

struct I::EchoCancellerMomentum
{
	Complex2D Input;
	Complex2D Loopback;
};

class EchoCancellerMomentum : public BaseFrequencyDomain<EchoCancellerMomentum, I::EchoCancellerMomentum>
{
	friend BaseFrequencyDomain<EchoCancellerMomentum, I::EchoCancellerMomentum>;
public:
	EchoSuppressionCovariance EchoSuppression;

private:
	struct Coefficients
	{
		float FilterbankRate = 125.f;
		float FilterLengthMilliseconds = 100.f;
		int NBands = 257;
		int NChannels = 2;
		int NChannelsLoopback = 2;
		int nBuffersInAnalysisFrame = 4;
		int nBuffersInSynthesisFrame = 4;
	} C;

	struct Parameters
	{
		ConstrainedType<float> LoopbackVarianceTConstant = { .180f, 0.001f, 60.f }; // seconds
		ConstrainedType<float> NearendVarianceTConstant = { .180f, 0.001f, 60.f }; // seconds
		ConstrainedType<float> CoefficientsVarianceTConstant = { .180f, 0.001f, 60.f }; // seconds
		ConstrainedType<float> NearendLimitdB = { -20.f, -100.f, 0.f }; // dB
		ConstrainedType<float> ChannelCorrelation = { 0.5f, 0.f, 1.f };
		bool FreezeUpdate = false;
	} P;

	struct Data
	{
		std::vector<std::vector<Eigen::ArrayXXcf>> Filters;
		int BufferLength;
		float LoopbackVarianceLambda, NearendVarianceLambda, CoefficientVarianceLambda, NearendLimit, ChannelCorrelation, ChannelCorrelation2;
		std::vector<std::vector<Eigen::ArrayXf>> CoefficientVariances, Momentums;
		std::vector<Eigen::ArrayXXcf> BuffersLoopback;
		std::vector<Eigen::ArrayXf> LoopbackVariances;
		Eigen::ArrayXcf LoopbackCovariance;
		Eigen::ArrayXf CorrelationCompensation;
		Eigen::ArrayXXf NearendVariance;
		int CircCounter;

		void Reset() 
		{
			for (auto &bufferLoopback : BuffersLoopback) { bufferLoopback.setZero(); }
			for (auto &filters : Filters) { for (auto& filter : filters) { filter.setZero(); } }
			for (auto &loopbackVariance : LoopbackVariances) { loopbackVariance.setConstant(100.f); }
			for (auto &coefficientVariances : CoefficientVariances) { for (auto& coefficientVariance : coefficientVariances) { coefficientVariance.setZero(); } }
			for (auto &momentums : Momentums) { for (auto& momentum : momentums) { momentum.setOnes(); } }
			CorrelationCompensation.setZero();
			LoopbackCovariance.setZero();
			NearendVariance.setConstant(10.f);
			CircCounter = 0;
		}
		bool InitializeMemory(const Coefficients& c) 
		{
			BufferLength = std::max(static_cast<int>(c.FilterLengthMilliseconds * c.FilterbankRate * 1e-3f) - (c.nBuffersInAnalysisFrame + c.nBuffersInSynthesisFrame - 1), 0)  + 1;
			BuffersLoopback.resize(c.NChannelsLoopback);
			for (auto &bufferLoopback : BuffersLoopback) { bufferLoopback.resize(c.NBands, BufferLength); }
			Filters.resize(c.NChannels);
			for (auto &filters : Filters) { filters.resize(c.NChannelsLoopback); for (auto& filter : filters) { filter.resize(c.NBands, BufferLength); } }
			LoopbackVariances.resize(c.NChannelsLoopback);
			for (auto &loopbackVariance : LoopbackVariances) { loopbackVariance.resize(c.NBands); }
			CoefficientVariances.resize(c.NChannels);
			for (auto& coefficientVariances : CoefficientVariances) { coefficientVariances.resize(c.NChannelsLoopback); for (auto& coefficientVariance : coefficientVariances) { coefficientVariance.resize(c.NBands); } }
			Momentums.resize(c.NChannels);
			for (auto& momentums : Momentums) { momentums.resize(c.NChannelsLoopback); for (auto& momentum : momentums) { momentum.resize(c.NBands); } }
			CorrelationCompensation.resize(c.NBands);
			LoopbackCovariance.resize(c.NBands);
			NearendVariance.resize(c.NBands, c.NChannels);
			return true;
		}
		size_t GetAllocatedMemorySize() const 
		{
			size_t size = 0;
			for (auto& filters : Filters) { for (auto& filter : filters) { size += filter.GetAllocatedMemorySize(); } }
			for (auto& buffer : BuffersLoopback) { size += buffer.GetAllocatedMemorySize(); }
			for (auto& loopback : LoopbackVariances) { size += loopback.GetAllocatedMemorySize(); }
			for (auto& coefs : CoefficientVariances) { for (auto& coef : coefs) { size += coef.GetAllocatedMemorySize(); } }
			for (auto& momentums : Momentums) { for (auto& momentum : momentums) { size += momentum.GetAllocatedMemorySize(); } }
			size += LoopbackCovariance.GetAllocatedMemorySize();
			size += NearendVariance.GetAllocatedMemorySize();
			size += CorrelationCompensation.GetAllocatedMemorySize();
			return size;
		}
		void OnParameterChange(const Parameters& p, const Coefficients& c) 
		{
			LoopbackVarianceLambda = 1.f - expf(-1.f / (c.FilterbankRate * p.LoopbackVarianceTConstant));
			NearendVarianceLambda = 1.f - expf(-1.f / (c.FilterbankRate * p.NearendVarianceTConstant));
			CoefficientVarianceLambda = 1.f - expf(-1.f / (c.FilterbankRate * p.CoefficientsVarianceTConstant));
			NearendLimit = powf(10.f, p.NearendLimitdB * 0.1f);
			ChannelCorrelation = p.ChannelCorrelation;
			ChannelCorrelation2 = ChannelCorrelation * ChannelCorrelation;
		}
	} D;

	DEFINEMEMBERALGORITHMS(1, EchoSuppression);

	auto InitializeMembers()
	{
		auto cES = EchoSuppression.GetCoefficients();
		cES.FilterbankRate = C.FilterbankRate;
		cES.NBands = C.NBands;
		cES.NChannels = C.NChannels;
		cES.NChannelsLoopback = C.NChannelsLoopback;
		return EchoSuppression.Initialize(cES);
	}

	void ProcessOn(Input xFreq, Output yFreq)
	{
		D.CircCounter = D.CircCounter <= 0 ? D.BufferLength - 1 : D.CircCounter - 1;
		// Put loopback signals into buffers
		for (auto i = 0; i < C.NChannelsLoopback; i++)
		{
			D.BuffersLoopback[i].col(D.CircCounter) = xFreq.Loopback.col(i); // put loopback signal i into buffer
		}
		//  filter and subtract from input
		Eigen::ArrayXXcf micEstimated = Eigen::ArrayXXcf::Zero(C.NBands, C.NChannels); // sum of loopback signals convolved with transfer function
		for (auto channel = 0; channel < C.NChannels; channel++)
		{
			for (auto chLoopback = 0; chLoopback < C.NChannelsLoopback; chLoopback++)
			{
				// estimated transferfunction x loopback
				for (auto j = 0; j < D.BufferLength - D.CircCounter; j++)
				{
					micEstimated.col(channel) += D.BuffersLoopback[chLoopback].col(j + D.CircCounter) * D.Filters[channel][chLoopback].col(j);
				}
				for (auto j = D.BufferLength - D.CircCounter; j < D.BufferLength; j++)
				{
					micEstimated.col(channel) += D.BuffersLoopback[chLoopback].col(j - D.BufferLength + D.CircCounter) * D.Filters[channel][chLoopback].col(j);
				}
			}
		}
		yFreq = xFreq.Input - micEstimated;
		D.NearendVariance += D.NearendVarianceLambda * (yFreq.abs2().max(micEstimated.abs2() * D.NearendLimit) - D.NearendVariance);

		if (!P.FreezeUpdate)
		{
			// update loopback (co-)variances
			for (auto chLoopback = 0; chLoopback < C.NChannelsLoopback; chLoopback++)
			{
				D.LoopbackVariances[chLoopback] += D.LoopbackVarianceLambda * (D.BuffersLoopback[chLoopback].col(D.CircCounter).abs2() - D.LoopbackVariances[chLoopback]); // variance of loopback signal i
			}

			if (C.NChannelsLoopback > 0)
			{
				D.LoopbackCovariance += D.LoopbackVarianceLambda * (D.BuffersLoopback[0].col(D.CircCounter).conjugate() * D.BuffersLoopback[C.NChannelsLoopback - 1].col(D.CircCounter) - D.LoopbackCovariance); // Covariance update
				D.CorrelationCompensation = 1.f - D.ChannelCorrelation2 * D.LoopbackCovariance.abs2() / (D.LoopbackVariances[0] * D.LoopbackVariances[C.NChannelsLoopback - 1] + 1e-20f);
			}

			for (auto channel = 0; channel < C.NChannels; channel++)
			{
				// update filters
				for (auto chLoopback = 0; chLoopback < C.NChannelsLoopback; chLoopback++)
				{
					Eigen::ArrayXcf signal1(C.NBands), signal2(C.NBands); // temporary
					signal1.setZero();
					signal2.setZero();
					signal2.real() = D.Momentums[channel][chLoopback] + D.BufferLength * D.CoefficientVariances[channel][chLoopback];
					D.Momentums[channel][chLoopback] = D.LoopbackVariances[chLoopback] * D.CorrelationCompensation;
					signal2.imag() = signal2.real() / (D.BufferLength * D.NearendVariance.col(channel) + (D.BufferLength + 2)*signal2.real()*D.Momentums[channel][chLoopback] + 1e-30f);
					signal2.imag() = signal2.imag().min(1.f / (D.Momentums[channel][chLoopback] * D.BufferLength + 1e-30f));
					D.Momentums[channel][chLoopback] = (1.f - signal2.imag() * D.Momentums[channel][chLoopback]) * signal2.real();
					signal1 = signal2.imag() * yFreq.col(channel);
					signal2 = D.ChannelCorrelation * D.LoopbackCovariance / (D.LoopbackVariances[C.NChannelsLoopback - 1 - chLoopback] + 1e-20f);
					D.Filters[channel][chLoopback].leftCols(D.BufferLength - D.CircCounter) += signal1.replicate(1, D.BufferLength - D.CircCounter)*(D.BuffersLoopback[chLoopback].rightCols(D.BufferLength - D.CircCounter).conjugate() - signal2.replicate(1, D.BufferLength - D.CircCounter)*D.BuffersLoopback[C.NChannelsLoopback - 1 - chLoopback].rightCols(D.BufferLength - D.CircCounter).conjugate());
					D.Filters[channel][chLoopback].rightCols(D.CircCounter) += signal1.replicate(1, D.CircCounter)*(D.BuffersLoopback[chLoopback].leftCols(D.CircCounter).conjugate() - signal2.replicate(1, D.CircCounter)*D.BuffersLoopback[C.NChannelsLoopback - 1 - chLoopback].leftCols(D.CircCounter).conjugate());

					D.CoefficientVariances[channel][chLoopback] += D.CoefficientVarianceLambda * ((signal1*(D.BuffersLoopback[chLoopback].col(D.CircCounter).conjugate() - signal2 * D.BuffersLoopback[C.NChannelsLoopback - 1 - chLoopback].col(D.CircCounter).conjugate())).abs2() - D.CoefficientVariances[channel][chLoopback]);
				}
			}
		}
		EchoSuppression.Process({ yFreq, xFreq.Loopback, micEstimated }, yFreq);
	}

	void ProcessOff(Input xFreq, Output yFreq) { yFreq = xFreq.Input; }
};