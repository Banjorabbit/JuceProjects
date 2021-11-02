#pragma once
#include "../BaseClasses/PureCRTP.h"
#include "GainCalculation.h"

struct I::EchoSuppressionCovariance
{
	Complex2D Input;
	Complex2D Loopback;
	Complex2D MicEstimated;
};

class EchoSuppressionCovariance : public BaseFrequencyDomain<EchoSuppressionCovariance, I::EchoSuppressionCovariance>
{
	friend BaseFrequencyDomain<EchoSuppressionCovariance, I::EchoSuppressionCovariance>;
public:
	EchoSuppressionCovariance() = default;
	VectorAlgo<GainCalculationSimple> GainCalculation;

private:
	struct Coefficients
	{
		float FilterbankRate = 125.f;
		int NBands = 257;
		int NChannels = 2;
		int NChannelsLoopback = 2;
	} C;

	struct Parameters
	{
		ConstrainedType<float> CovarianceTConstant = { 5.f, 0.001f, 60.f }; // seconds
		ConstrainedType<float> LoopbackPowerTConstant = { 0.150f, 0.001f, 60.f }; // seconds
		ConstrainedType<float> SnrSmoother = { 0.98f, 0.f, 1.f };
		ConstrainedType<float> SpeechDistortionRatio = { 1.f, 0.f, 100.f };
		ConstrainedType<float> GainMinimumdB = { -30.f, -100.f, 0.f }; // dB
		ConstrainedType<float> GainAttackTConstant = { .001f, 1e-6f, 60.f }; // seconds
		ConstrainedType<float> GainReleaseTConstant = { 0.02f, 1e-6f, 60.f }; // seconds
	} P;

	struct Data
	{
		Eigen::ArrayXf LoopbackPower, LoopbackVariance;
		Eigen::ArrayXXcf Covariance;
		float LoopbackPowerLambda, CovarianceLambda;
		float GainMinimum, GainAttackLambda, GainReleaseLambda;

		void Reset() 
		{
			Covariance.setConstant(1e-10f);
			LoopbackVariance.setConstant(1e-10f);
			LoopbackPower.setZero();
		}
		bool InitializeMemory(const Coefficients& c) 
		{
			LoopbackVariance.resize(c.NBands);
			LoopbackPower.resize(c.NBands);
			Covariance.resize(c.NBands, c.NChannels);
			return true;
		}
		size_t GetAllocatedMemorySize() const 
		{ 
			size_t size = 0;
			size += LoopbackPower.GetAllocatedMemorySize(); 
			size += LoopbackVariance.GetAllocatedMemorySize();
			size += Covariance.GetAllocatedMemorySize();
			return size;
		}
		void OnParameterChange(const Parameters& p, const Coefficients& c) 
		{
			LoopbackPowerLambda = 1.f - expf(-1.f / (c.FilterbankRate * p.LoopbackPowerTConstant));
			CovarianceLambda = 1.f - expf(-1.f / (c.FilterbankRate * p.CovarianceTConstant));
			GainMinimum = powf(10.f, p.GainMinimumdB * .05f);
			GainAttackLambda = 1.f - expf(-1.f / (c.FilterbankRate*p.GainAttackTConstant));
			GainReleaseLambda = 1.f - expf(-1.f / (c.FilterbankRate*p.GainReleaseTConstant));
		}
	} D;

	DEFINEMEMBERALGORITHMS(1, GainCalculation)

	auto InitializeMembers() 
	{
		GainCalculation.resize(C.NChannels);
		auto flag = true;
		for (auto& gainCalculation : GainCalculation)
		{
			auto setup = gainCalculation.GetSetup();
			setup.Coefficients.FilterbankRate = C.FilterbankRate;
			setup.Coefficients.NBands = C.NBands;
			setup.Coefficients.NChannels = 1;
			setup.Parameters.Method = setup.Parameters.Smoothing;
			setup.Parameters.MinimumdB = -30;
			flag &= gainCalculation.Initialize(setup);
		}
		return flag;
	}

	void ProcessOn(Input xFreq, Output yFreq)
	{
		Eigen::ArrayXf power = (xFreq.Loopback.matrix() * Eigen::VectorXf::Ones(C.NChannelsLoopback)).array().abs2(); // multiplying with VectorXf:Ones is faster than rowwise sum
		D.LoopbackPower += D.LoopbackPowerLambda * (power - D.LoopbackPower);
		D.LoopbackVariance += (power > 1e-10f).select(D.CovarianceLambda * (power - D.LoopbackVariance), 0.f);
		Eigen::ArrayXf snr = D.LoopbackVariance / power.max(D.LoopbackPower).max(1e-30f);
		Eigen::ArrayXXf gain(C.NBands, C.NChannels);
		for (auto channel = 0; channel < xFreq.Input.cols(); channel++)
		{
			D.Covariance.col(channel) += (power > 1e-10f).select(D.CovarianceLambda * (xFreq.Input.col(channel)*xFreq.MicEstimated.col(channel).conjugate() - D.Covariance.col(channel)), 0.f);
			Eigen::ArrayXf snrAposteriori = xFreq.Input.col(channel).abs2() * snr / D.Covariance.col(channel).abs();
			GainCalculation[channel].Process(snrAposteriori, gain.col(channel));
		}
		yFreq = xFreq.Input * gain;
	}

	void ProcessOff(Input xFreq, Output yFreq) { yFreq = xFreq.Input; }
};