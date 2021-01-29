#pragma once
#include "../BaseClasses/PureCRTP.h"
#include "../FrequencyDomain/DetectTonal.h"
#include "../Filterbank.h"
#include "../FFT.h"

class SeparateTonal : public AsynchronousBase<SeparateTonal>
{
	friend Base<SeparateTonal>;

public:
	DetectTonal TonalDetector;
	FilterbankAnalysis Filterbank;
	FilterbankSynthesis FilterbankInverse;
	FFTReal FFT;
	FFTReal FFTHalf;


	int GetLatencySamples() const { return D.LatencySamples; }

	auto GetAsynchronousCoefficients(const std::vector<int> bufferSizesSuggested) const
	{
		Coefficients c = C;
		auto bufferSize = static_cast<int>(powf(2, std::round(std::log2(.01f*c.SampleRate))));
		auto fftSize = 16 * bufferSize;
		if (FFTReal::IsFFTSizeValid(fftSize))
		{
			c.BufferSize = bufferSize;
			c.FrameSize = 8 * bufferSize;
		}
		return c;
	}

private:
	struct Coefficients
	{
		float SampleRate = 44100.f;
		int BufferSize = 512;
		int FrameSize = 4096;
		int NChannelsIn = 2;
		AsynchronousBufferType AsynchronousBuffer = VARIABLE_SIZE;
	} C;
	struct Parameters {} P;

	struct Data
	{
		int LatencySamples, FFTSize, NBands, NBandsHalf;
		void Reset() { }
		bool InitializeMemory(const Coefficients& c)
		{
			LatencySamples = c.FrameSize - c.BufferSize;
			FFTSize = 2 * c.FrameSize;
			NBands = FFTSize / 2 + 1;
			NBandsHalf = FFTSize / 4 + 1;
			return true;
		}
		size_t GetAllocatedMemorySize() const { return 0; }
		void OnParameterChange(const Parameters& p, const Coefficients& c) { }
	} D;

	DEFINEMEMBERALGORITHMS(5, TonalDetector, Filterbank, FilterbankInverse, FFT, FFTHalf)

	auto InitializeMembers()
	{
		auto sT = TonalDetector.GetSetup();
		sT.Coefficients.FilterbankRate = C.SampleRate / C.BufferSize;
		sT.Coefficients.NBands = D.NBandsHalf;
		sT.Coefficients.NChannels = C.NChannelsIn;
		sT.Coefficients.SampleRate = C.SampleRate;
		sT.Coefficients.WindowSizeFreqHz = 200;
		sT.Parameters.PowerCompression = 0.3f;
		sT.Parameters.TonalTConstant = 0.15f;
		sT.Parameters.TonalThreshold = 2.f;
		sT.Parameters.TransientTConstant = 0.005f;
		sT.Parameters.TransientThreshold = 1.25f;
		auto flag = TonalDetector.Initialize(sT);

		auto sF = Filterbank.GetSetup();
		sF.Coefficients.NChannels = C.NChannelsIn;
		sF.Coefficients.FrameSize = C.FrameSize;
		sF.Coefficients.FFTSize = D.FFTSize;
		sF.Coefficients.BufferSize = C.BufferSize;
		sF.Parameters.WindowType = sF.Parameters.HannWindow;
		flag &= Filterbank.Initialize(sF);

		auto sFI = FilterbankInverse.GetSetup();
		sFI.Coefficients.NChannels = C.NChannelsIn;
		sFI.Coefficients.FrameSize = D.FFTSize;
		sFI.Coefficients.FFTSize = D.FFTSize;
		sFI.Coefficients.BufferSize = C.BufferSize;
		sFI.Parameters.WindowType = sFI.Parameters.Rectangular;
		sFI.Parameters.Gain = static_cast<float>(2 * C.BufferSize) / C.FrameSize;
		flag &= FilterbankInverse.Initialize(sFI);

		auto cFFT = FFT.GetCoefficients();
		cFFT.FFTSize = D.FFTSize;
		flag &= FFT.Initialize(cFFT);

		auto cFFTHalf = FFTHalf.GetCoefficients();
		cFFTHalf.FFTSize = D.FFTSize / 2;
		flag &= FFTHalf.Initialize(cFFTHalf);

		return flag;
	}

	void CalculateGain(const Eigen::Ref<const Eigen::ArrayXXf>& input, Eigen::Ref<Eigen::ArrayXXcf> output)
	{
		Eigen::ArrayXXf h(D.FFTSize, C.NChannelsIn);
		FFTHalf.Inverse(input.cast<std::complex<float>>(), h.topRows(D.FFTSize/2));
		h.bottomRows(D.FFTSize * 3 / 4).setZero();
		FFT.Process(h, output);
	}

	void ProcessOn(Input xTime, Output yTime)
	{
		Eigen::ArrayXXcf xFreq(D.NBands, C.NChannelsIn);
		Filterbank.Process(xTime, xFreq);
		Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> DecisionTonal(D.NBandsHalf, C.NChannelsIn);
		Eigen::ArrayXXcf xFreqHalf = Eigen::Map<Eigen::ArrayXXcf, 0, Eigen::InnerStride<2> >(xFreq.data(), D.NBandsHalf, C.NChannelsIn);
		TonalDetector.Process(xFreqHalf, DecisionTonal);
		Eigen::ArrayXXcf Filter(D.NBands, C.NChannelsIn);
		CalculateGain(DecisionTonal.cast<float>(), Filter);
		xFreq *= Filter;
		FilterbankInverse.Process(xFreq, yTime);
	}

	void ProcessOff(Input xTime, Output yTime) { yTime = xTime; }
};