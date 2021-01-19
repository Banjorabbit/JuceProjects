#pragma once
#include "../BaseClasses/PureCRTP.h"

struct I::InterpolateTonal
{
	RealX2 xEnergy; // 2 energy spectrums to interpolate between.
	RealX2 xPhase; // 2 phase vectors corresponding to xEnergy. This format is chosen so xPhase.col(0) doesn't have to be recomputed if it was xPhase.col(1) last time.
	float PointFractional; // floats in range [0,1[, which is a fractional interpolation point
};

class InterpolateTonal : public BaseFrequencyDomain<InterpolateTonal, I::InterpolateTonal, O::Complex>
{
	friend BaseFrequencyDomain<InterpolateTonal, I::InterpolateTonal, O::Complex>;

	struct Coefficients 
	{
		int NBands = 1023;
	} C;
	struct Parameters {} P;

	struct Data 
	{
		Eigen::ArrayXf Phase;
		void Reset()
		{
			Phase.setZero();
		}
		bool InitializeMemory(const Coefficients& c)
		{
			Phase.resize(c.NBands);
			return true;
		}
		size_t GetAllocatedMemorySize() const
		{
			return Phase.GetAllocatedMemorySize();
		}
		void OnParameterChange(const Parameters& p, const Coefficients& c) { }
	} D;

	void ProcessOn(Input xFreq, Output yFreq) 
	{
		Eigen::ArrayXf energy = (xFreq.xEnergy.col(0) + xFreq.PointFractional * (xFreq.xEnergy.col(1) - xFreq.xEnergy.col(0)));
		D.Phase += xFreq.xPhase.col(1) - xFreq.xPhase.col(0);

		Eigen::Array<bool, Eigen::Dynamic, 1> localMax(C.NBands);
		localMax(0) = false;
		localMax(C.NBands - 1) = false;
		localMax.segment(1, C.NBands - 2) = energy.head(C.NBands - 2) < energy.segment(1, C.NBands - 2) && energy.tail(C.NBands - 2) < energy.segment(1, C.NBands - 2); // size NBands - 2

		// detect tones leaked to neighbour bins and match phase accross frequency
		auto lastEnergy = 0.f;
		for (auto i = 1; i < C.NBands - 1; i++)
		{
			if (localMax(i))
			{
				D.Phase(i + 1) = D.Phase(i) + static_cast<const float>(PI);
				if (energy(i) > lastEnergy) { D.Phase(i - 1) = D.Phase(i) - static_cast<const float>(PI); }
				lastEnergy = energy(i);
			}
		}
		for (auto i = 1; i < C.NBands - 1; i++)
		{
			// wrap to +- PI
			if (D.Phase(i) > static_cast<const float>(PI)) { D.Phase(i) -= 2.f*static_cast<const float>(PI); }
			else if (D.Phase(i) < -static_cast<const float>(PI)) { D.Phase(i) += 2.f*static_cast<const float>(PI); }
		}

		// with angle
		energy = energy.sqrt();
		yFreq.real() = energy * D.Phase.cos();
		yFreq.imag() = energy * D.Phase.sin();
	}

	void ProcessOff(Input xFreq, Output yFreq) { yFreq.setZero(); }
};