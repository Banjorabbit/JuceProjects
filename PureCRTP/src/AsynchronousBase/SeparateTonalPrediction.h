#pragma once
#include "../BaseClasses/PureCRTP.h"
#include "../Filterbank.h"
#include "../FrequencyDomain/EchoCancellerNLMSMinPow.h"

class SeparateTonalPrediction : public AsynchronousBase<SeparateTonalPrediction>
{
	friend Base<SeparateTonalPrediction>;

public:
	FilterbankAnalysis Filterbank;
	FilterbankSynthesis FilterbankInverse;
	VectorAlgo<EchoCancellerNLMSMinPow> Predictor; // nChannels in EchoCancellerNLMS assumes they share RX channel, but since they are different here, we need a vector of predictors
	
	int GetLatencySamples() const { return D.LatencySamples; }

	auto GetAsynchronousCoefficients(const std::vector<int> bufferSizesSuggested) const
	{
		Coefficients c = C;
		auto bufferSize = static_cast<int>(powf(2, std::round(std::log2(.01f*c.SampleRate))));
		auto fftSize = 4 * bufferSize;
		if (FFTReal::IsFFTSizeValid(fftSize))
		{
			c.BufferSize = bufferSize;
			c.FrameSize = fftSize;
		}
		return c;
	}

private:
	struct Coefficients
	{
		float SampleRate = 44100.f;
		int BufferSize = 512;
		int FrameSize = 2048;
		int NChannelsIn = 2;
		int PredictionDelay = 1;
		AsynchronousBufferType AsynchronousBuffer = VARIABLE_SIZE;
	} C;

	struct Parameters { } P;

	struct Data
	{
		Eigen::ArrayXXcf xFreqDelayed;
		int LatencySamples, NBands;
		void Reset() { xFreqDelayed.setZero(); }
		bool InitializeMemory(const Coefficients& c)
		{
			LatencySamples = c.FrameSize - c.BufferSize;
			NBands = c.FrameSize / 2 + 1;
			xFreqDelayed.resize(NBands, c.NChannelsIn);
			return true;
		}
		size_t GetAllocatedMemorySize() const { return xFreqDelayed.GetAllocatedMemorySize(); }
		void OnParameterChange(const Parameters& p, const Coefficients& c) { }
	} D;

	DEFINEMEMBERALGORITHMS(3, Filterbank, FilterbankInverse, Predictor)

	auto InitializeMembers()
	{
		auto FilterbankRate = C.SampleRate / C.BufferSize;

		auto sF = Filterbank.GetSetup();
		sF.Coefficients.NChannels = C.NChannelsIn;
		sF.Coefficients.FrameSize = C.FrameSize;
		sF.Coefficients.FFTSize = C.FrameSize;
		sF.Coefficients.BufferSize = C.BufferSize;
		sF.Parameters.WindowType = sF.Parameters.HannWindow;
		auto flag = Filterbank.Initialize(sF);

		auto sFI = FilterbankInverse.GetSetup();
		sFI.Coefficients.NChannels = C.NChannelsIn;
		sFI.Coefficients.FrameSize = C.FrameSize;
		sFI.Coefficients.FFTSize = C.FrameSize;
		sFI.Coefficients.BufferSize = C.BufferSize;
		sFI.Parameters.WindowType = sFI.Parameters.HannWindow;
		sFI.Parameters.Gain = 2.f/3.f;
		flag &= FilterbankInverse.Initialize(sFI);

		Predictor.resize(C.NChannelsIn);
		auto sP = Predictor[0].GetSetup(); // all channels are setup the same
		sP.Coefficients.FilterbankRate = FilterbankRate;
		sP.Coefficients.NFilters = 5;
		sP.Coefficients.FilterLength[0] = 1;
		sP.Coefficients.FilterLength[1] = 2;
		sP.Coefficients.FilterLength[2] = 4;
		sP.Coefficients.FilterLength[3] = 8;
		sP.Coefficients.FilterLength[4] = 16;
		sP.Coefficients.NBands = D.NBands;
		sP.Coefficients.NChannels = 1; // nChannels is encoded in the size of std::vector<EchoCancellerNLMS>
		sP.Parameters.NearendLimitdB = -40;
		flag &= Predictor.Initialize(sP);

		// calculate synthesis gain
		Eigen::ArrayXf ASWindow = Filterbank.GetWindow() * FilterbankInverse.GetWindow();
		float gain = Eigen::Map<Eigen::ArrayXXf>(ASWindow.data(), C.BufferSize, C.FrameSize / C.BufferSize).rowwise().sum().mean();
		sFI.Parameters.Gain = 1.f / gain;
		FilterbankInverse.SetParameters(sFI.Parameters);

		return flag;
	}

	void ProcessOn(Input xTime, Output yTime)
	{
		Eigen::ArrayXXcf xFreq(D.NBands, C.NChannelsIn), yFreq(D.NBands, C.NChannelsIn);

		Filterbank.Process(xTime, xFreq);
		for (auto channel = 0; channel < C.NChannelsIn; channel++)
		{
			Predictor[channel].Process({ xFreq.col(channel), D.xFreqDelayed.col(channel) }, yFreq.col(channel));
		}
		D.xFreqDelayed = xFreq;
		FilterbankInverse.Process(xFreq - yFreq, yTime);
	}

	void ProcessOff(Input xTime, Output yTime) { yTime = xTime; }
};