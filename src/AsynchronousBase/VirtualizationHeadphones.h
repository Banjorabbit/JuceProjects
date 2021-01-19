#pragma once
#include "../BaseClasses/PureCRTP.h"
#include "FFTConvolution.h"
#include "../Upsampling2XCubic.h"

class VirtualizationHeadphones : public AsynchronousBase<VirtualizationHeadphones>
{
	friend Base<VirtualizationHeadphones>;

public:
	int GetNChannelsOut() const { return 2; }
	int GetLatencySamples() const { return D.DelaySamples; }

	auto GetAsynchronousCoefficients(const std::vector<int> bufferSizesSuggested) const
	{
		Coefficients c = C;

		int scaleFactor = 1;
		if (c.SampleRate == c.Hz88200 || c.SampleRate == c.Hz96000) { scaleFactor = 2; }
		else if (c.SampleRate == c.Hz176400 || c.SampleRate == c.Hz192000) { scaleFactor = 4; }

		auto bufferSizes = bufferSizesSuggested;

		// vector of bufferSizes to try
		bufferSizes.push_back(std::max(static_cast<int>(std::pow(2, std::round(std::log2(bufferSizes[0] + FilterSize))) - FilterSize), 1)); // round to nearest int of log2(size + Algo.FilterSize)
		bufferSizes.push_back(std::max(static_cast<int>(std::pow(2, std::ceil(std::log2(bufferSizes[0] + FilterSize))) - FilterSize), 1)); // round to ceiling int of log2(size + Algo.FilterSize)
		bufferSizes.push_back(256);
		bufferSizes.push_back(32);

		// return with first value that is a valid size
		for (auto& size : bufferSizes)
		{
			if (FFTReal::IsFFTSizeValid(size + FilterSize * scaleFactor))
			{
				c.BufferSize = size;
				return c;
			}
		}
		return c; // failed to find valid BufferSize so use default
	}

	static constexpr int FilterSize = 256; // FilterSize is hardcoded to 256 since only those coefficients exist
	VectorAlgo<FFTConvolution> ConvolveFilter; // vector of FFConvolutions. One for each input channel
	UpsamplingPower2Cubic Upsampling; // this is used to upsample loaded filters

private:
	struct Coefficients
	{
		ConstrainedType<int> NChannelsIn = { 2, 1, 6 }; // if EnabledSub == true, then the last channel will be the Sub channel
		int BufferSize = 256;
		enum SampleRates { Hz44100, Hz48000, Hz88200, Hz96000, Hz176400, Hz192000};
		SampleRates SampleRate = Hz44100;
		enum VirtualChannels { FrontLeft, FrontRight, BackLeft, BackRight, Center };
		VirtualChannels DirectionsChannels[5] = { FrontLeft, FrontRight, BackLeft, BackRight, Center}; // Only the first C.NChannels DirectionsChannels are used (if EnabledSub == true, then C.NChannels-1 are used)
		bool EnabledSub = false;
		AsynchronousBufferType AsynchronousBuffer = VARIABLE_SIZE;
	} C;

	struct Parameters {} P;

	struct Data 
	{
		int ScaleFactor;
		int DelaySamples;
		int NChannelsVirtual;
		Eigen::ArrayXf BufferSub;
		void Reset() { BufferSub.setZero(); }
		bool InitializeMemory(const Coefficients& c) 
		{ 
			if (c.EnabledSub) { NChannelsVirtual = c.NChannelsIn - 1; }
			else { NChannelsVirtual = std::min(5, static_cast<int>(c.NChannelsIn)); } // only the first 5 channels have virtual directions (size of C.DirectionsChannels)

			ScaleFactor = 1;
			if (c.SampleRate == c.Hz176400 || c.SampleRate == c.Hz192000) { ScaleFactor = 4; }
			else if (c.SampleRate == c.Hz88200 || c.SampleRate == c.Hz96000) { ScaleFactor = 2; }
			DelaySamples = 31 * ScaleFactor; // hardcoded based on impulse responses in .cpp file
			BufferSub.resize(DelaySamples);
			return true; 
		}
		size_t GetAllocatedMemorySize() const { return BufferSub.GetAllocatedMemorySize(); }
		void OnParameterChange(const Parameters& p, const Coefficients& c) { }
	} D;

	void LoadFilter(Coefficients::VirtualChannels channel, Eigen::Ref<Eigen::ArrayXXf> filter); // defined in .cpp file

	DEFINEMEMBERALGORITHMS(2, ConvolveFilter, Upsampling);

	auto InitializeMembers()
	{
		auto filterSizeScaled = FilterSize * D.ScaleFactor;

		ConvolveFilter.resize(D.NChannelsVirtual);
		auto setup = ConvolveFilter[0].GetSetup();
		setup.Coefficients.BufferSize = C.BufferSize;
		setup.Coefficients.FilterSize = filterSizeScaled;
		setup.Coefficients.NChannelsIn = 1;
		setup.Coefficients.NFiltersPerChannel = 2;
		auto flag = ConvolveFilter.Initialize(setup);
		flag &= Upsampling.Initialize();
		if (flag)
		{
			// calculate filter in ConvolveFIlter
			for (auto i = 0;i < ConvolveFilter.size();i++)
			{
				Eigen::ArrayXXf filter(filterSizeScaled, 2);
				LoadFilter(C.DirectionsChannels[i], filter);
				ConvolveFilter[i].SetPersistentInput(filter);
			}
		}
		return flag;
	}

	void ProcessOn(Input xTime, Output yTime) 
	{
		Eigen::ArrayXXf output(C.BufferSize, 2);
		ConvolveFilter[0].Process( xTime.col(0), yTime);
		for (auto i = 1;i < ConvolveFilter.size();i++)
		{
			ConvolveFilter[i].Process( xTime.col(i), output);
			yTime += output;
		}
		if (C.EnabledSub) 
		{ 
			yTime.topRows(D.DelaySamples) += D.BufferSub.replicate(1, 2);
			yTime.bottomRows(C.BufferSize - D.DelaySamples) += xTime.rightCols<1>().head(C.BufferSize - D.DelaySamples).replicate(1,2);
			D.BufferSub = xTime.rightCols<1>().tail(D.DelaySamples);
		}
	}

	void ProcessOff(Input xTime, Output yTime) { yTime = xTime.leftCols<1>().replicate(1,2); }
};
