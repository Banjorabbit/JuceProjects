#pragma once
#include "BaseClasses/PureCRTP.h"

// Estimates the spectral mass using the zero-crossing rate in the time-domain
//
// author: Kristian Timm Andersen

class SpectralMass : public Base<SpectralMass>
{
	friend Base<SpectralMass>;

	struct Coefficients 
	{
		float SampleRate = 44100.f;
		int NChannels = 2;
	} C;

	struct Parameters 
	{
		float SmoothTConstant = 0.01f;
	} P;

	struct Data 
	{
		Eigen::ArrayXf XOld, SpectralMassHz, ZeroCrossOld, SpectralMassOld;
		float SmoothLambda;
		void Reset() 
		{
			XOld.setZero();
			SpectralMassOld.setConstant(1000.f);
			SpectralMassHz.setConstant(1000.f);
			ZeroCrossOld.setZero();
		}
		bool InitializeMemory(const Coefficients& c)
		{
			XOld.resize(c.NChannels);
			SpectralMassOld.resize(c.NChannels);
			SpectralMassHz.resize(c.NChannels);
			ZeroCrossOld.resize(c.NChannels);
			return true;
		}
		size_t GetAllocatedMemorySize() const 
		{
			size_t size = SpectralMassHz.GetAllocatedMemorySize();
			size += ZeroCrossOld.GetAllocatedMemorySize();
			size += SpectralMassOld.GetAllocatedMemorySize();
			size += XOld.GetAllocatedMemorySize();
			return size; 
		}
		void OnParameterChange(const Parameters& p, const Coefficients& c) 
		{
			SmoothLambda = 1.f - expf(-1.f / (c.SampleRate * p.SmoothTConstant));
		}
	} D;

	void ProcessOn(Input xTime, Output spectralMassHz) 
	{
		for (auto channel = 0; channel < C.NChannels; channel++)
		{
			for (auto i = 0; i < xTime.rows(); i++)
			{
				// detect distance between each zero-crossing
				if (xTime(i, channel) * D.XOld(channel) < 0.f) // if zero-crossing
				{
					const float p = D.XOld(channel) / (D.XOld(channel) - xTime(i, channel) + 1e-20f) + i - 1.f; // linear interpolation between values to determine zero crossing place
					const float ZeroCrossDistance = p - D.ZeroCrossOld(channel);
					D.ZeroCrossOld(channel) = p;
					D.SpectralMassHz(channel) = std::min(C.SampleRate / (2.f * ZeroCrossDistance + 1e-20f), C.SampleRate / 2.f); // instantanious spectral mass in Hz (2 zero-crossings is 1 period)
				}
				
				const float weightedLambda = D.SmoothLambda * std::fabs(xTime(i, channel)); // weight by abs(x) so more smoothing is done when silence
				D.SpectralMassOld(channel) += weightedLambda * (D.SpectralMassHz(channel) - D.SpectralMassOld(channel)); // 1st-order low pass filter
				spectralMassHz(i, channel) = D.SpectralMassOld(channel);

				D.XOld(channel) = xTime(i, channel);
			}
			D.ZeroCrossOld(channel) -= xTime.rows(); // update so next ProcessOn() call has correct value
		}
	}

	void ProcessOff(Input xTime, Output spectralMassHz) { spectralMassHz.setZero(); }
};