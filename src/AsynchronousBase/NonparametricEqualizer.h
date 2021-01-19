#pragma once
#include "../BaseClasses/PureCRTP.h"
#include "../FFT.h"
#include "../DesignFIRNonParametric.h"

// Nonparametric equalizer
//
// Given a set of frequencies and corresponding gains, a filter is calculated and applied to the input signal.
//
// The implementation is hardcoded to use Hann windows with 50% overlap between the frames.
//
// author: Kristian Timm Andersen

struct I::NonparametricEqualizerPersistent
{
	Real Frequencies; // set of frequencies
	Real GaindB; // set of corresponding gains in dB
};

class NonparametricEqualizer : public AsynchronousBase<NonparametricEqualizer, I::Real2D, O::Real2D, I::NonparametricEqualizerPersistent>
{
	friend Base<NonparametricEqualizer, I::Real2D, O::Real2D, I::NonparametricEqualizerPersistent>;

public:
	FFTReal FFT;
	DesignFIRNonParametric FilterCalculation;

	auto GetAsynchronousCoefficients(const std::vector<int> bufferSizesSuggested) const
	{
		auto c = GetCoefficients();
		auto bufferSizes = bufferSizesSuggested;
		bufferSizes.push_back(std::max(static_cast<int>(std::pow(2, std::round(std::log2(bufferSizesSuggested[0] + c.FilterSize))) - c.FilterSize), 1)); // round to nearest int of log2(size + c.Filter)
		bufferSizes.push_back(std::max(static_cast<int>(std::pow(2, std::ceil(std::log2(bufferSizesSuggested[0] + c.FilterSize))) - c.FilterSize), 1)); // round to ceiling int of log2(size + c.Filter)
		bufferSizes.push_back(32);

		// return with first value that is a valid size
		for (auto& size : bufferSizes)
		{
			if (FFTReal::IsFFTSizeValid(size + c.FilterSize))
			{
				c.BufferSize = size;
				return c; // return with first match
			}
		}
		return c; // failed to find valid BufferSize so use default
	}

private:
	struct Coefficients 
	{
		int NChannelsIn = 2;
		int BufferSize = 256;
		int FilterSize = 384;
		float SampleRate = 48e3;
		AsynchronousBufferType AsynchronousBuffer = VARIABLE_SIZE;
	} C;

	struct Parameters {} P;

	struct Data
	{
		Eigen::ArrayXXf TimeBuffer;
		Eigen::ArrayXcf FrequencyFilter;
		int FFTSize, NBands;
		void Reset() 
		{
			TimeBuffer.setZero();
			FrequencyFilter.setZero();
		}
		bool InitializeMemory(const Coefficients& c) 
		{
			FFTSize = c.BufferSize + c.FilterSize;
			NBands = FFTSize / 2 + 1;
			TimeBuffer.resize(FFTSize, c.NChannelsIn);
			FrequencyFilter.resize(NBands);
			return true;
		}
		size_t GetAllocatedMemorySize() const 
		{
			size_t size = TimeBuffer.GetAllocatedMemorySize();
			size += FrequencyFilter.GetAllocatedMemorySize();
			return size;
		}
		void OnParameterChange(const Parameters& p, const Coefficients& c)	{}
	} D;

	DEFINEMEMBERALGORITHMS(2, FFT, FilterCalculation);

	auto InitializeMembers()
	{
		auto sFFT = FFT.GetSetup();
		sFFT.Coefficients.FFTSize = D.FFTSize;
		auto flag = FFT.Initialize(sFFT);
		auto sFC = FilterCalculation.GetSetup();
		sFC.Coefficients.FilterSize = C.FilterSize;
		sFC.Coefficients.SampleRate = C.SampleRate;
		flag &= FilterCalculation.Initialize(sFC);
		return flag;
	}

	void ProcessPersistentInput(InputPersistent input)
	{
		Eigen::ArrayXf timeBuffer(D.FFTSize);
		FilterCalculation.Process({ input.Frequencies, input.GaindB }, timeBuffer.head(C.FilterSize));
		timeBuffer.tail(D.FFTSize - C.FilterSize).setZero();
		FFT.Process(timeBuffer, D.FrequencyFilter);
	}

	void ProcessOn(Input xTime, Output yTime) 
	{
		Eigen::ArrayXf timeBuffer(D.FFTSize);
		for (auto channel = 0; channel < C.NChannelsIn; channel++)
		{
			Eigen::ArrayXcf frequencyBuffer(D.NBands);
			timeBuffer.head(C.BufferSize) = xTime.col(channel);
			timeBuffer.tail(C.FilterSize).setZero();
			FFT.Process(timeBuffer, frequencyBuffer);
			frequencyBuffer *= D.FrequencyFilter;
			FFT.Inverse(frequencyBuffer, timeBuffer);
			D.TimeBuffer.col(channel) += timeBuffer;
			yTime.col(channel) = D.TimeBuffer.block(0, channel, C.BufferSize, 1);
			D.TimeBuffer.block(0, channel, C.FilterSize, 1) = D.TimeBuffer.block(C.BufferSize, channel, C.FilterSize, 1);
			D.TimeBuffer.block(C.FilterSize, channel, C.BufferSize, 1).setZero();
		}	
	}
	
	void ProcessOff(Input xTime, Output yTime) { yTime = xTime; }
};