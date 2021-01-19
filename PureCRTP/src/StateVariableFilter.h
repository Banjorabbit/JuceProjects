#pragma once
#include "BaseClasses/PureCRTP.h"

struct O::StateVariableFilter
{
	Real2D HighPass;
	Real2D LowPass;
	Real2D BandPass;
};

class StateVariableFilter : public Base<StateVariableFilter, I::Real2D, O::StateVariableFilter>
{
	friend Base<StateVariableFilter, I::Real2D, O::StateVariableFilter>;

	struct Coefficients
	{
		float SampleRate = 16e3f;
		int NChannels = 2;
		int BufferSize = 128;
	} C;

	struct Parameters
	{
		ConstrainedType<float> Cutoff = { 1000.f, 1.f, 20e3f };
		ConstrainedType<float> Q = { static_cast<float>(BUTTERWORTH_Q), 1e-3f, 100.f };
	} P;

	struct Data
	{
		Eigen::ArrayXf State1, State2;
		float c0, c1, c2, c3;
		void Reset()
		{
			State1.setZero();
			State2.setZero();
		}
		bool InitializeMemory(const Coefficients& c)
		{
			State1.resize(c.NChannels);
			State2.resize(c.NChannels);
			return true;
		}
		size_t GetAllocatedMemorySize() const
		{
			auto size = State1.GetAllocatedMemorySize();
			size += State2.GetAllocatedMemorySize();
			return size;
		}
		void OnParameterChange(const Parameters& p, const Coefficients& c)
		{
			const auto g = static_cast<float>(std::tan(PI*p.Cutoff / c.SampleRate));
			c1 = g / p.Q;
			c2 = g * p.Q;
			c3 = c2 + 1.f;
			c0 = 1.f / (1.f + c1 * (c2 + 1.f));
		}
	} D;

	void ProcessOn(Input xTime, Output yTime)
	{
		for (auto sample = 0; sample < C.BufferSize; sample++)
		{
			for (auto channel = 0; channel < C.NChannels; channel++)
			{
				yTime.HighPass(sample, channel) = D.c0 * (xTime(sample, channel) - D.State2(channel) - D.c3 * D.State1(channel));
				const auto x1 = D.c1 * yTime.HighPass(sample, channel);
				yTime.BandPass(sample, channel) = x1 + D.State1(channel);
				const auto x2 = D.c2 * yTime.BandPass(sample, channel);
				yTime.LowPass(sample, channel) = x2 + D.State2(channel);

				D.State1(channel) = x1 + yTime.BandPass(sample, channel);
				D.State2(channel) = x2 + yTime.LowPass(sample, channel);
			}

		}
		// combine to make Notch, Bell, HighShelf, LowShelf:
		// D.g0 = 10.^(gain./20);
		// D.g1 = D.g0 - 1;
		//yTime.Notch = xTime - yTime.BandPass;
		//yTime.Bell = xTime + D.g1 * yTime.BandPass;
		//yTime.HighShelf = D.g0 * (yTime.BandPass + yTime.HighPass) + yTime.LowPass;
		//yTime.LowShelf = D.g0 * (yTime.BandPass + yTime.LowPass) + yTime.HighPass;
	}

	void ProcessOff(Input xTime, Output yTime)
	{
		yTime.HighPass = xTime;
		yTime.LowPass = xTime;
		yTime.BandPass = xTime;
	}
};