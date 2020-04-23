#pragma once
#pragma warning (disable : 4127 ) //4127 removes "unused variable" warnings from Eigen
#define EIGEN_MPL2_ONLY // don't allow LGPL licensed code from Eigen
#include <Eigen/Dense> // Eigen Library. Potentially redundant since this is also included in InputOutput.h
#pragma warning (default : 4127) // set "unused variable" warning to default warning level
#include "InputOutput.h" 

// AsynchronousStreaming is a wrapper using CRTP for an algorithm that derives from Base. 
// It allows streaming data asynchronously in and out of the algorithm, which is particularly 
// useful for processing audio data from a host that sends/receives audio frames that are either 
// not supported by the algorithm or produces sub-optimal results.
//
// author: Kristian Timm Andersen

// Talgo is assumed to derive from Base<Talgo>
template<typename Tderived, typename Talgo>
class AsynchronousStreaming
{
public:
	Talgo Algo;

	size_t GetAllocatedMemorySize() const
	{
		return BufferIn.GetAllocatedMemorySize() + BufferOut.GetAllocatedMemorySize() + Algo.GetAllocatedMemorySize();
	}

	static size_t GetStaticMemorySize() { return sizeof(Tderived) + Talgo::GetStaticMemorySize(); }

	int GetLatencyTotalSamples() const { return LatencyTotalSamples; }
	bool GetSynchronousProcessing() const { return SynchronousProcessing; }
	int GetBufferSizeInternal() const { return BufferSizeInternal; }
	int GetNChannelsOut() const { return NChannelsOut; }
	int GetNChannels() const { return NChannels; }
	int GetUpdatedCoefficients(decltype(Algo.GetCoefficients())& c, const int bufferSizeExpected, const int nChannels, const float sampleRate) const
	{
		std::vector<int> bufferSizesSuggested = { bufferSizeExpected, std::max(bufferSizeExpected / 2, 1), std::max(bufferSizeExpected / 4, 1), std::max(bufferSizeExpected / 8, 1),
										 std::max(static_cast<int>(std::pow(2, std::round(std::log2(bufferSizeExpected)))), 1),
										 std::max(static_cast<int>(std::pow(2, std::ceil(std::log2(bufferSizeExpected)))), 1) };
		return static_cast<const Tderived&>(*this).UpdateCoefficients(c, sampleRate, nChannels, bufferSizeExpected, bufferSizesSuggested);
	}

	void Reset()
	{
		ResetData();
		Algo.Reset();
	}

	bool Initialize(const int bufferSizeExpected, const int nChannels, const float sampleRate)
	{
		return Initialize(bufferSizeExpected, nChannels, sampleRate, Algo.GetCoefficients());
	}

	bool Initialize(const int bufferSizeExpected, const int nChannels, const float sampleRate, const decltype(Algo.GetCoefficients()) coefficients)
	{
		auto& derived = static_cast<Tderived&>(*this);

		Algo.SetCoefficients(coefficients);

		auto c = Algo.GetCoefficients();
		BufferSizeInternal = GetUpdatedCoefficients(c, bufferSizeExpected, nChannels, sampleRate);
		if (BufferSizeInternal <= 0) { return false; }
		bool successFlag = Algo.Initialize(c);

		NChannels = nChannels;
		NChannelsOut = derived.CalculateNChannelsOut(NChannels); // some algorithms does not care how many channels it gets and therefore does know how many input channels it is set up to process

		if (std::round(static_cast<float>(bufferSizeExpected) / BufferSizeInternal) == static_cast<float>(bufferSizeExpected) / BufferSizeInternal)
		{
			SynchronousProcessing = true;
		}
		else { SynchronousProcessing = false; }

		BufferIn.resize(BufferSizeInternal, NChannels);
		BufferOut.resize(BufferSizeInternal, NChannelsOut);
		ResetData();

		if (SynchronousProcessing) { LatencyTotalSamples = derived.CalculateLatencySamples(); }
		else { LatencyTotalSamples = BufferSizeInternal + derived.CalculateLatencySamples(); }

		return successFlag;
	}

	void Process(typename Talgo::Input xTime, typename Talgo::Output yTime)
	{
		if (SynchronousProcessing) { ProcessEntireBuffers(xTime, yTime); }
		else { ProcessSmallerBuffers(xTime, yTime); }
	}
private:

	// Process the entire xTime buffer and put result in yTime. This function does not result in a delay.
	void ProcessEntireBuffers(typename Talgo::Input xTime, typename Talgo::Output yTime)
	{
		int i = 0;
		// Process as many full buffers as possible
		for (i; i <= xTime.rows() - BufferSizeInternal; i += BufferSizeInternal)
		{
			Algo.Process(xTime.middleRows(i, BufferSizeInternal), yTime.middleRows(i, BufferSizeInternal));
		}
		// if we have been given a size that is not an integer multiple of BufferSizeInternal, zeropad and process. This is used if host does not give the expected length.
		int remainingSamples = std::min(static_cast<int>(xTime.rows()) - i, BufferSizeInternal);
		for (i; i < xTime.rows(); i += remainingSamples)
		{
			BufferIn.topRows(remainingSamples) = xTime.middleRows(i, remainingSamples);
			BufferIn.bottomRows(BufferSizeInternal - remainingSamples).setZero();
			Algo.Process(BufferIn, BufferOut);
			yTime.middleRows(i, remainingSamples) = BufferOut.topRows(remainingSamples);
			remainingSamples = std::min(static_cast<int>(xTime.rows()) - i, BufferSizeInternal);
		}
	}

	// Copy xTime buffer into BufferIn and process when it is full. While copying from xTime, also copy BufferOut (result of last Process) to yTime. This results in a delay equal to length of BufferSizeInternal.
	// This function is used when bufferSizeExpected is not an integer multiple of BufferSizeInternal (for instance if bufferSizeExpected is smaller than BufferSizeInternal, or 1.5x, 2.5x,... times larger).
	void ProcessSmallerBuffers(typename Talgo::Input xTime, typename Talgo::Output yTime)
	{
		for (auto i = 0; i < xTime.rows(); i++)
		{
			BufferIn.row(Index) = xTime.row(i);
			yTime.row(i) = BufferOut.row(Index);
			Index++;
			if (Index == BufferSizeInternal)
			{
				Algo.Process(BufferIn, BufferOut);
				Index = 0;
			}
		}
	}

	void ResetData()
	{
		Index = 0;
		BufferIn.setZero();
		BufferOut.setZero();
		BufferSizeEmpty = BufferSizeInternal;
	}

	// this is set in Initialize() 
	int BufferSizeInternal = 0;
	bool SynchronousProcessing = true;
	int Index = 0;
	int NChannels = 0;
	int NChannelsOut = 0;
	int LatencyTotalSamples = 0;
	int BufferSizeEmpty = 0;
	Eigen::ArrayXXf BufferIn, BufferOut;
};