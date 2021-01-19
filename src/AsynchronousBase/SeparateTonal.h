#pragma once
#include "../BaseClasses/PureCRTP.h"
#include "../FrequencyDomain/DetectTonal.h"
#include "../Filterbank.h"

class SeparateTonal : public AsynchronousBase<SeparateTonal>
{
	friend Base<SeparateTonal>;

public:
	DetectTonal TonalDetector;
	FilterbankAnalysis Filterbank;
	FilterbankSynthesis FilterbankInverse;

	int GetLatencySamples() const { return D.LatencySamples; }

	auto GetAsynchronousCoefficients(const std::vector<int> bufferSizesSuggested) const
	{
		Coefficients c = C;
		auto bufferSize = static_cast<int>(powf(2, std::round(std::log2(.01f*c.SampleRate))));
		auto fftSize = 8 * bufferSize;
		if (FFTReal::IsFFTSizeValid(fftSize))
		{
			c.BufferSize = bufferSize;
			c.FFTSize = fftSize;
		}
		return c;
	}

private:
	struct Coefficients
	{
		float SampleRate = 44100.f;
		int BufferSize = 512;
		int FFTSize = 4096;
		int NChannelsIn = 2;
		AsynchronousBufferType AsynchronousBuffer = VARIABLE_SIZE;
	} C;
	struct Parameters {} P;

	struct Data
	{
		int LatencySamples, NBands;
		void Reset() { }
		bool InitializeMemory(const Coefficients& c)
		{
			LatencySamples = c.FFTSize - c.BufferSize;
			NBands = c.FFTSize / 2 + 1;
			return true;
		}
		size_t GetAllocatedMemorySize() const { return 0; }
		void OnParameterChange(const Parameters& p, const Coefficients& c) { }
	} D;

	DEFINEMEMBERALGORITHMS(3, TonalDetector, Filterbank, FilterbankInverse)

	auto InitializeMembers()
	{
		auto sT = TonalDetector.GetSetup();
		sT.Coefficients.FilterbankRate = C.SampleRate / C.BufferSize;
		sT.Coefficients.NBands = D.NBands;
		sT.Coefficients.NChannels = C.NChannelsIn;
		sT.Coefficients.SampleRate = C.SampleRate;
		sT.Coefficients.WindowSizeFreqHz = 200;
		sT.Parameters.PowerCompression = 0.3f;
		sT.Parameters.TonalTConstant = 0.15f;
		sT.Parameters.TonalThreshold = 2.f;
		auto flag = TonalDetector.Initialize(sT);

		auto sF = Filterbank.GetSetup();
		sF.Coefficients.NChannels = C.NChannelsIn;
		sF.Coefficients.FrameSize = C.FFTSize;
		sF.Coefficients.FFTSize = C.FFTSize;
		sF.Coefficients.BufferSize = C.BufferSize;
		sF.Parameters.WindowType = sF.Parameters.HannWindow;
		flag &= Filterbank.Initialize(sF);

		auto sFI = FilterbankInverse.GetSetup();
		sFI.Coefficients.NChannels = C.NChannelsIn;
		sFI.Coefficients.FrameSize = C.FFTSize;
		sFI.Coefficients.FFTSize = C.FFTSize;
		sFI.Coefficients.BufferSize = C.BufferSize;
		sFI.Parameters.WindowType = sFI.Parameters.HannWindow;
		sFI.Parameters.Gain = static_cast<float>(4 * C.BufferSize) / C.FFTSize;
		flag &= FilterbankInverse.Initialize(sFI);

		return flag;
	}

	void ProcessOn(Input xTime, Output yTime)
	{
		Eigen::ArrayXXcf xFreq(D.NBands, C.NChannelsIn);
		Filterbank.Process(xTime, xFreq);
		Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> DecisionTonal(D.NBands, C.NChannelsIn);
		TonalDetector.Process(xFreq, DecisionTonal);
		xFreq *= DecisionTonal.cast<float>();
		FilterbankInverse.Process(xFreq, yTime);
	}

	void ProcessOff(Input xTime, Output yTime) { yTime = xTime; }
};