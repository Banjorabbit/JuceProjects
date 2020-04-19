#pragma once
#include "../BaseClasses/PureCRTP.h"
#include "FFTConvolution.h"
#include "../Upsampling2XCubic.h"

class VirtualizationHeadphones : public Base<VirtualizationHeadphones>
{
	friend Base<VirtualizationHeadphones>;

public:
	static constexpr int FilterSize = 256; // FilterSize is hardcoded to 256 since only those coefficients exist
	int GetLatencySamples() const { return D.DelaySamples; }
	VectorAlgo<FFTConvolution> ConvolveFilter; // vector of FFConvolutions. One for each input channel
	UpsamplingPower2Cubic Upsampling; // this is used to upsample loaded filters

private:
	struct Coefficients
	{
		ConstrainedType<int> NChannels = { 2, 1, 6 }; // if EnabledSub == true, then the last channel will be the Sub channel
		int BufferSize = 256;
		enum SampleRates { Hz44100, Hz48000, Hz88200, Hz96000, Hz176400, Hz192000};
		SampleRates SampleRate = Hz44100;
		enum VirtualChannels { FrontLeft, FrontRight, BackLeft, BackRight, Center };
		VirtualChannels DirectionsChannels[5] = { FrontLeft, FrontRight, BackLeft, BackRight, Center}; // Only the first C.NChannels DirectionsChannels are used (if EnabledSub == true, then C.NChannels-1 are used)
		bool EnabledSub = false;
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
			if (c.EnabledSub) { NChannelsVirtual = c.NChannels - 1; }
			else { NChannelsVirtual = std::min(5, static_cast<int>(c.NChannels)); } // only the first 5 channels have virtual directions (size of C.DirectionsChannels)

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
		setup.Coefficients.NChannels = 1;
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

class VirtualizationHeadphonesStreaming : public AsynchronousStreaming< VirtualizationHeadphonesStreaming, VirtualizationHeadphones>
{
	friend AsynchronousStreaming<VirtualizationHeadphonesStreaming, VirtualizationHeadphones>;

	int GetLatencySamples(const decltype(Algo.GetCoefficients())& c) const { return Algo.GetLatencySamples(); }
	int GetNChannelsOut(const int nChannels) const { return 2; }

	int UpdateCoefficients(decltype(Algo.GetCoefficients())& c, const float sampleRate, const int nChannels, const int bufferSizeExpected, std::vector<int> bufferSizesSuggested) const
	{
		int scaleFactor = 1;
		if (sampleRate == 44100.f) { c.SampleRate = c.Hz44100; }
		else if (sampleRate == 48000.f) { c.SampleRate = c.Hz48000; }
		else if (sampleRate == 88200.f) { c.SampleRate = c.Hz88200; scaleFactor = 2; }
		else if (sampleRate == 96000.f) { c.SampleRate = c.Hz96000; scaleFactor = 2; }
		else if (sampleRate == 176200.f) { c.SampleRate = c.Hz176400; scaleFactor = 4; }
		else if (sampleRate == 192000.f) { c.SampleRate = c.Hz192000; scaleFactor = 4; }
		else { return -1; }
		c.NChannels = nChannels;

		// vector of bufferSizes to try
		bufferSizesSuggested.push_back(std::max(static_cast<int>(std::pow(2, std::round(std::log2(bufferSizeExpected + Algo.FilterSize))) - Algo.FilterSize), 1)); // round to nearest int of log2(size + Algo.FilterSize)
		bufferSizesSuggested.push_back(std::max(static_cast<int>(std::pow(2, std::ceil(std::log2(bufferSizeExpected + Algo.FilterSize))) - Algo.FilterSize), 1)); // round to ceiling int of log2(size + Algo.FilterSize)
		bufferSizesSuggested.push_back(256);
		bufferSizesSuggested.push_back(32);
		bufferSizesSuggested.push_back(c.BufferSize);

		// return with first value that is a valid size
		for (auto& size : bufferSizesSuggested)
		{
			if (FFTReal::IsFFTSizeValid(size + Algo.FilterSize * scaleFactor))
			{
				c.BufferSize = size;
				return size;
			}
		}
		// failed to set BufferSize correctly
		return -1;
	}

	
};