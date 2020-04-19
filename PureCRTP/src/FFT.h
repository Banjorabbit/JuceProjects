#pragma once
#include "BaseClasses/PureCRTP.h"
#include "../Utilities/pffft.h"

// Wrapper for real pffft. 
//
// pffft requires 16 byte alligned data, and we don't know the allocation of input/output (float is 4 bytes). 
// The assert in pffft.c is: #define VALIGNED(ptr) ((((uintptr_t)(ptr)) & 0xF) == 0)
// Especially, frequency domain data has length FFTSize/2+1 so if first channel is 16 byte alligned then second channel is guaranteed not to be 16 byte alligned.
// For this reason, FFT transform is using memory allocated on the stack.
//
// This class has a public inverse FFT function.
//
// author: Kristian Timm Andersen
class FFTReal : public Base<FFTReal, I::Real2D, O::Complex2D>
{
	friend Base<FFTReal, I::Real2D, O::Complex2D>;

public:
	
	void Inverse(const Eigen::Ref<const Eigen::ArrayXXcf>& xFreq, Eigen::Ref<Eigen::ArrayXXf> yTime) const 
	{
		if (GetEnabled() && GetInitialized())
		{
			for (auto channel = 0; channel < xFreq.cols(); channel++)
			{
				yTime(0,channel) = xFreq(0, channel).real();
				yTime(1,channel) = xFreq(C.FFTSize / 2, channel).real();
				std::memcpy(&yTime(2,channel), xFreq.col(channel).data()+1, (C.FFTSize - 2) * sizeof(float));
				pffft_transform_ordered(D.Setup.get(), yTime.col(channel).data(), yTime.col(channel).data(), nullptr, PFFFT_BACKWARD);
				yTime.col(channel) *= D.Scale;
			}
		}
		else { yTime.setZero();	}
	}

	static bool IsFFTSizeValid(const int FFTSize)
	{
		if (FFTSize % 32 != 0 || FFTSize < 32) { return false; } // first check size ís integer factor of 32
		PFFFT_Setup *setup = pffft_new_setup(FFTSize, PFFFT_REAL);
		if (setup == nullptr) { return false; }
		pffft_destroy_setup(setup);
		return true;
	}

private:
	struct Coefficients
	{
		ConstrainedType<int> FFTSize = { 512, 32, 16384 }; // pffft requires minimum 32 real size
	} C;

	struct Parameters { } P;

	struct Data
	{
		float Scale;
		std::shared_ptr<PFFFT_Setup> Setup;
		void Reset() { }
		bool InitializeMemory(const Coefficients& c)
		{
			Setup = std::shared_ptr<PFFFT_Setup>(pffft_new_setup(c.FFTSize, PFFFT_REAL), pffft_destroy_setup);
			assert(Setup.get() != nullptr);
			if (Setup.get() == nullptr) { return false; }
			return true;
		}
		size_t GetAllocatedMemorySize() const
		{
			if (Setup.get() == nullptr) { return 0; }
			else { return pffft_get_setup_size(Setup.get()); }
		}
		void OnParameterChange(const Parameters& p, const Coefficients& c) { Scale = 1.f / static_cast<float>(c.FFTSize); }
	} D;

	void ProcessOn(Input xTime, Output yFreq) const 
	{
		// yFreq is not 16 byte alligned due to FFTSize/2+1 size so we can't write output to it from FFT transform
		for (auto channel = 0; channel < xTime.cols(); channel++)
		{
			Eigen::ArrayXf out(C.FFTSize);
			pffft_transform_ordered(D.Setup.get(), xTime.col(channel).data(), out.data(), nullptr, PFFFT_FORWARD);
			yFreq(0, channel) = out(0);
			yFreq(C.FFTSize / 2, channel) = out(1);
			std::memcpy(&yFreq.real()(1, channel), &out(2), (C.FFTSize - 2) * sizeof(float));
		}		
	}

	void ProcessOff(Input xTime, Output yFreq) { yFreq.setZero(); }
};

// Wrapper for inverse real pffft. 
//
// pffft requires 16 byte alligned data, and we don't know the allocation of input/output (float is 4 bytes). 
// The assert in pffft.c is: #define VALIGNED(ptr) ((((uintptr_t)(ptr)) & 0xF) == 0)
// Especially, frequency domain data has length FFTSize/2+1 so if first channel is 16 byte alligned then second channel is guaranteed not to be 16 byte alligned.
//
// author: Kristian Timm Andersen
class FFTRealInverse : public Base<FFTRealInverse, I::Complex2D, O::Real2D>
{
	friend Base<FFTRealInverse, I::Complex2D, O::Real2D>;

public:
	bool IsFFTSizeValid(const int FFTSize)
	{
		PFFFT_Setup *setup = pffft_new_setup(FFTSize, PFFFT_REAL);
		if (setup == nullptr) { return false; }
		pffft_destroy_setup(setup);
		return true;
	}

private:
	struct Coefficients
	{
		ConstrainedType<int> FFTSize = { 512, 32, 16384 }; // pffft requires minimum 32 real size
	} C;

	struct Parameters { } P;

	struct Data
	{
		float Scale;
		std::shared_ptr<PFFFT_Setup> Setup;
		void Reset() { }
		bool InitializeMemory(const Coefficients& c)
		{
			Setup = std::shared_ptr<PFFFT_Setup>(pffft_new_setup(c.FFTSize, PFFFT_REAL), pffft_destroy_setup);
			assert(Setup.get() != nullptr);
			if (Setup.get() == nullptr) { return false; }
			return true;
		}
		size_t GetAllocatedMemorySize() const
		{
			if (Setup.get() == nullptr) { return 0; }
			else { return pffft_get_setup_size(Setup.get()); }
		}
		void OnParameterChange(const Parameters& p, const Coefficients& c) { Scale = 1.f / static_cast<float>(c.FFTSize); }
	} D;

	void ProcessOn(Input xFreq, Output yTime) const
	{
		for (auto channel = 0; channel < xFreq.cols(); channel++)
		{
			yTime(0, channel) = xFreq(0, channel).real();
			yTime(1, channel) = xFreq(C.FFTSize / 2, channel).real();
			std::memcpy(&yTime(2, channel), xFreq.col(channel).data() + 1, (C.FFTSize - 2) * sizeof(float));
			pffft_transform_ordered(D.Setup.get(), yTime.col(channel).data(), yTime.col(channel).data(), nullptr, PFFFT_BACKWARD);
			yTime.col(channel) *= D.Scale;
		}
	}

	void ProcessOff(Input xFreq, Output yTime) { yTime.setZero(); }
};