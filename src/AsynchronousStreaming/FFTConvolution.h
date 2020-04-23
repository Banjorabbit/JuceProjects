#pragma once
#include "../BaseClasses/PureCRTP.h"
#include "../Filterbank.h"
#include "../FFT.h"

class FFTConvolution : public Base<FFTConvolution>
{
	friend Base<FFTConvolution>;

public:
	FilterbankAnalysis Filterbank;
	VectorAlgo<FilterbankSynthesis> FilterbankInverse;
	FFTReal FFT;

private:

	struct Coefficients 
	{
		int BufferSize = 256;
		int NChannels = 8;
		int NFiltersPerChannel = 1;
		int FilterSize = 512;
	} C;

	struct Parameters { } P;

	struct Data {
		int FFTSize, NBands;
		Eigen::ArrayXXcf FilterFreq;
		void Reset() { }
		bool InitializeMemory(const Coefficients& c)	
		{
			FFTSize = c.BufferSize + c.FilterSize;
			NBands = FFTSize / 2 + 1;
			FilterFreq.resize(NBands, c.NFiltersPerChannel);
			FilterFreq.setOnes();
			return true;
		}
		size_t GetAllocatedMemorySize() const { return FilterFreq.GetAllocatedMemorySize(); }
		void OnParameterChange(Parameters& p, const Coefficients& c) { }
	} D;

	DEFINEMEMBERALGORITHMS(3, Filterbank, FilterbankInverse, FFT);

	auto InitializeMembers()
	{
		auto s = Filterbank.GetSetup();
		s.Coefficients.BufferSize = C.BufferSize;
		s.Coefficients.FFTSize = D.FFTSize;
		s.Coefficients.FrameSize = C.BufferSize;
		s.Coefficients.NChannels = C.NChannels;
		s.Parameters.WindowType = s.Parameters.Rectangular;
		bool flag = Filterbank.Initialize(s);

		FilterbankInverse.resize(C.NFiltersPerChannel);
		auto sInverse = FilterbankInverse[0].GetSetup();
		sInverse.Coefficients.BufferSize = C.BufferSize;
		sInverse.Coefficients.FFTSize = D.FFTSize;
		sInverse.Coefficients.FrameSize = D.FFTSize;
		sInverse.Coefficients.NChannels = C.NChannels;
		sInverse.Parameters.WindowType = sInverse.Parameters.Rectangular;
		flag &= FilterbankInverse.Initialize(sInverse);

		auto cFFT = FFT.GetCoefficients();
		cFFT.FFTSize = D.FFTSize;
		flag &= FFT.Initialize(cFFT);
		return flag;
	}

	void ProcessPersistentInput(InputPersistent filter) 
	{
		// update filter in frequency domain
		Eigen::ArrayXXf timeFilter(D.FFTSize, C.NFiltersPerChannel);
		timeFilter.topRows(C.FilterSize) = filter;
		timeFilter.bottomRows(D.FFTSize - C.FilterSize).setZero();
		FFT.Process(timeFilter, D.FilterFreq);
	}

	void ProcessOn(Input xTime, Output yTime) 
	{
		// FFT convolution
		Eigen::ArrayXXcf xFreq(D.NBands, C.NChannels);
		Filterbank.Process(xTime, xFreq);
		for (auto filter = 0;filter < C.NFiltersPerChannel;filter++)
		{
			Eigen::ArrayXXcf yFreq(D.NBands, C.NChannels);
			yFreq = xFreq * D.FilterFreq.col(filter).replicate(1, C.NChannels);
			FilterbankInverse[filter].Process(yFreq, yTime.middleCols(filter * C.NChannels, C.NChannels));
		}
	}

	void ProcessOff(Input xTime, Output yTime) { yTime = xTime; }
};

class FFTConvolutionStreaming : public AsynchronousStreaming<FFTConvolutionStreaming, FFTConvolution>
{
	friend AsynchronousStreaming<FFTConvolutionStreaming, FFTConvolution>;

	int CalculateLatencySamples() const { return 0; }
	int CalculateNChannelsOut(const int nChannels) const { const auto c = Algo.GetCoefficients(); return c.NChannels * c.NFiltersPerChannel; }

	int UpdateCoefficients(decltype(Algo.GetCoefficients())& c, const float sampleRate, const int nChannels, const int bufferSizeExpected, std::vector<int> bufferSizesSuggested) const
	{
		c.NChannels = nChannels;
		// vector of bufferSizes to try
		bufferSizesSuggested.push_back(std::max(static_cast<int>(std::pow(2, std::round(std::log2(bufferSizeExpected + c.FilterSize))) - c.FilterSize), 1)); // round to nearest int of log2(size + c.Filter)
		bufferSizesSuggested.push_back(std::max(static_cast<int>(std::pow(2, std::ceil(std::log2(bufferSizeExpected + c.FilterSize))) - c.FilterSize), 1)); // round to ceiling int of log2(size + c.Filter)
		bufferSizesSuggested.push_back(32);
		bufferSizesSuggested.push_back(c.BufferSize);

		// return with first value that is a valid size
		for (auto& size : bufferSizesSuggested)
		{
			if (FFTReal::IsFFTSizeValid(size + c.FilterSize))
			{
				c.BufferSize = size;
				return c.BufferSize;
			}
		}
		// failed to find valid BufferSize
		return -1;
	}
};