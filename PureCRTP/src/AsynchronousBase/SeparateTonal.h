#pragma once
#include "../BaseClasses/PureCRTP.h"
#include "../FrequencyDomain/DetectTonalPowerMin.h"
#include "../FrequencyDomain//DetectTransientSmooth.h"
#include "../Filterbank.h"
#include "../FrequencyDomain/EchoCancellerNLMS.h"
#include "../CircBuffer.h"

class SeparateTonal : public AsynchronousBase<SeparateTonal>
{
	friend Base<SeparateTonal>;

public:
	DetectTonalPowerMin TonalDetector;
	DetectTransientSmooth TransientDetector;
	FilterbankAnalysis Filterbank;
	FilterbankAnalysis FilterbankTransient;
	FilterbankAnalysis FilterbankDelayed;
	FilterbankSynthesis FilterbankInverse;
	VectorAlgo<EchoCancellerNLMS> Predictor; // nChannels in EchoCancellerNLMS assumes they share RX channel, but since they are different here, we need a vector of predictors
	CircBuffer Delay;


	int GetLatencySamples() const { return D.LatencySamples; }

	auto GetAsynchronousCoefficients(const std::vector<int> bufferSizesSuggested) const
	{
		Coefficients c = C;
		auto bufferSize = static_cast<int>(powf(2, std::round(std::log2(.01f*c.SampleRate))));
		auto fftSize = 8 * bufferSize;
		auto fftSizeTransient = 2 * bufferSize;
		if (FFTReal::IsFFTSizeValid(fftSize) & FFTReal::IsFFTSizeValid(fftSizeTransient))
		{
			c.BufferSize = bufferSize;
			c.FrameSize = fftSize;
			c.FrameSizeTransient = fftSizeTransient;
		}
		return c;
	}

private:
	struct Coefficients
	{
		float SampleRate = 44100.f;
		int BufferSize = 512;
		int FrameSize = 4096;
		int FrameSizeTransient = 1024;
		int NChannelsIn = 2;
		int PredictionDelay = 1280;
		AsynchronousBufferType AsynchronousBuffer = VARIABLE_SIZE;
	} C;

	struct Parameters { } P;

	struct Data
	{
		int LatencySamples, NBands;
		void Reset() { }
		bool InitializeMemory(const Coefficients& c)
		{
			LatencySamples = c.FrameSize - c.BufferSize;
			NBands = c.FrameSize / 2 + 1;
			return true;
		}
		size_t GetAllocatedMemorySize() const { return 0; }
		void OnParameterChange(const Parameters& p, const Coefficients& c) { }
	} D;

	DEFINEMEMBERALGORITHMS(8, TonalDetector, TransientDetector, Filterbank, FilterbankTransient, FilterbankDelayed, FilterbankInverse, Predictor, Delay)

	auto InitializeMembers()
	{
		auto FilterbankRate = C.SampleRate / C.BufferSize;

		auto sT = TonalDetector.GetSetup();
		sT.Coefficients.FilterbankRate = FilterbankRate;
		sT.Coefficients.NBands = D.NBands;
		sT.Coefficients.NChannels = C.NChannelsIn;
		sT.Coefficients.SampleRate = C.SampleRate;
		sT.Coefficients.WindowSizeFreqHz = 200;
		sT.Parameters.PowerCompression = 0.3f;
		sT.Parameters.TonalTConstant = 0.1f;
		sT.Parameters.TonalThreshold = 2.f;
		auto flag = TonalDetector.Initialize(sT);

		auto sTT = TransientDetector.GetSetup();
		sTT.Coefficients.FilterbankRate = FilterbankRate;
		sTT.Coefficients.NBands = D.NBands;
		sTT.Coefficients.NChannels = C.NChannelsIn;
		sTT.Coefficients.SampleRate = C.SampleRate;
		sTT.Coefficients.WindowSizeFreqHz = 800.f;
		sTT.Parameters.TConstant = 0.1f;
		sTT.Parameters.Threshold = 1.5f;
		sTT.Parameters.PowerCompression = 0.3f;
		flag &= TransientDetector.Initialize(sTT);

		auto sF = Filterbank.GetSetup();
		sF.Coefficients.NChannels = C.NChannelsIn;
		sF.Coefficients.FrameSize = C.FrameSize;
		sF.Coefficients.FFTSize = C.FrameSize;
		sF.Coefficients.BufferSize = C.BufferSize;
		sF.Parameters.WindowType = sF.Parameters.HannWindow;
		flag &= Filterbank.Initialize(sF);

		auto sFT = FilterbankTransient.GetSetup();
		sFT.Coefficients.NChannels = C.NChannelsIn;
		sFT.Coefficients.FrameSize = C.FrameSize; // keep using FrameSize but set a shorter window
		sFT.Coefficients.FFTSize = C.FrameSize;
		sFT.Coefficients.BufferSize = C.BufferSize;
		sFT.Parameters.WindowType = sFT.Parameters.UserDefined; 
		sFT.Parameters.Gain = 4; // window is 4 times short than Filterbank, so gain should be 4 times higher
		flag &= FilterbankTransient.Initialize(sFT);
		// set userdefined window
		Eigen::ArrayXf windowTransient(C.FrameSize);
		int wStart = (C.FrameSize - C.FrameSizeTransient) / 2;
		windowTransient.head(wStart).setZero();
		windowTransient.segment(wStart, C.FrameSizeTransient) = FilterbankAnalysis::GetHannWindow(C.FrameSizeTransient);
		windowTransient.tail(wStart).setZero();
		FilterbankTransient.SetWindow(windowTransient);
		
		auto sFD = FilterbankDelayed.GetSetup();
		sFD.Coefficients.NChannels = C.NChannelsIn;
		sFD.Coefficients.FrameSize = C.FrameSize;
		sFD.Coefficients.FFTSize = C.FrameSize;
		sFD.Coefficients.BufferSize = C.BufferSize;
		sFD.Parameters.WindowType = sFD.Parameters.HannWindow;
		flag &= FilterbankDelayed.Initialize(sFD);

		auto sFI = FilterbankInverse.GetSetup();
		sFI.Coefficients.NChannels = C.NChannelsIn;
		sFI.Coefficients.FrameSize = C.FrameSize;
		sFI.Coefficients.FFTSize = C.FrameSize;
		sFI.Coefficients.BufferSize = C.BufferSize;
		sFI.Parameters.WindowType = sFI.Parameters.HannWindow;
		sFI.Parameters.Gain = 1;
		flag &= FilterbankInverse.Initialize(sFI);

		Predictor.resize(C.NChannelsIn);
		auto sP = Predictor[0].GetSetup(); // all channels are setup the same
		sP.Coefficients.FilterbankRate = FilterbankRate;
		sP.Coefficients.FilterLength = 1;
		sP.Coefficients.NBands = D.NBands;
		sP.Coefficients.NChannels = 1; // nChannels is encoded in the size of std::vector<EchoCancellerNLMS>
		sP.Parameters.NearendLimitdB = -40;
		flag &= Predictor.Initialize(sP);

		auto cD = Delay.GetCoefficients();
		cD.DelayLength = C.PredictionDelay;
		cD.NChannels = C.NChannelsIn;
		flag &= Delay.Initialize(cD);

		// calculate synthesis gain
		Eigen::ArrayXf ASWindow = Filterbank.GetWindow() * FilterbankInverse.GetWindow();
		float gain = Eigen::Map<Eigen::ArrayXXf>(ASWindow.data(), C.BufferSize, C.FrameSize / C.BufferSize).rowwise().sum().mean();
		sFI.Parameters.Gain = 1.f/gain;
		FilterbankInverse.SetParameters(sFI.Parameters);

		return flag;
	}

	void ProcessOn(Input xTime, Output yTime)
	{
		Eigen::ArrayXXf xDelay(C.BufferSize, C.NChannelsIn);
		Eigen::ArrayXXcf xFreq(D.NBands, C.NChannelsIn), xFreqDelayed(D.NBands, C.NChannelsIn), xFreqTransient(D.NBands, C.NChannelsIn);
		Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> DecisionTonal(D.NBands, C.NChannelsIn), DecisionTransient(D.NBands, C.NChannelsIn);

		Delay.Process(xTime, xDelay);
		FilterbankDelayed.Process(xDelay, xFreqDelayed);
		Filterbank.Process(xTime, xFreq);
		FilterbankTransient.Process(xTime, xFreqTransient);

		TransientDetector.Process(xFreqTransient, DecisionTransient);
		TonalDetector.Process(xFreq, DecisionTonal);

		
		Eigen::ArrayXXcf FiltersOptimal = xFreq / (xFreqDelayed + 1e-10f);
		Eigen::ArrayXXcf yFreq(D.NBands,C.NChannelsIn);
		for (auto channel = 0; channel < C.NChannelsIn; channel++)
		{
			std::vector<Eigen::ArrayXXcf> filters = Predictor[channel].GetFilters();
			for (auto iband = 0; iband < D.NBands; iband++)
			{
				if (DecisionTransient(iband, channel)) { filters[0](0, iband) = FiltersOptimal(iband, channel) * 1e-6f; }
				if (DecisionTonal(iband, channel)) { filters[0](0, iband) = FiltersOptimal(iband, channel); }
			}
			Predictor[channel].SetFilters(filters);

			Predictor[channel].Process({ xFreq.col(channel), xFreqDelayed.col(channel) }, yFreq.col(channel)); // predictor actually outputs the non-periodic part!
		}
		yFreq = (xFreq.abs2() < yFreq.abs2()).select(xFreq, yFreq); // only use predictor if energy is less or equal to input

		FilterbankInverse.Process(xFreq - yFreq, yTime);
	}

	void ProcessOff(Input xTime, Output yTime) { yTime = xTime; }
};