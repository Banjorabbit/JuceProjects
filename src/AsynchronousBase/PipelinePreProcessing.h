#pragma once
#include "../BaseClasses/PureCRTP.h"
#include "IIR2ndDF.h"
#include "../Filterbank.h"
#include "../FrequencyDomain/EchoCancellerMomentum.h"
#include "../FrequencyDomain/DetectVoiceActivation.h"
#include "../FrequencyDomain/BeamformerAdaptive.h"
#include "../FrequencyDomain/NoiseSuppression.h"

class PipelinePreProcessing : public AsynchronousBase<PipelinePreProcessing>
{
	friend Base<PipelinePreProcessing>;
	struct Coefficients; // forward declaration so CalculateNBuffersInFrame(Coefficients& c) can be defined

public:
	int GetNChannelsOut() const { return 1; }
	auto GetLatencySamples() const { return D.LatencySamples; }
	auto GetFFTSize() const { return D.FFTSize; }
	static int CalculateNBuffersInFrame(const Coefficients& c)
	{
		int nBuffersInFrame = 1 + static_cast<int>(c.DelayDesired * c.SampleRate / c.BufferSize);
		nBuffersInFrame = std::max(static_cast<int>(nBuffersInFrame / 4) * 4, 4); // force factor of 4
		return nBuffersInFrame;
	}

	auto GetAsynchronousCoefficients(const std::vector<int> bufferSizesSuggested) const
	{
		Coefficients c = C;
		auto bufferSizes = bufferSizesSuggested;

		// start suggestion should be size 2^x close to 10ms
		int goodSize = static_cast<int>(powf(2, std::round(std::log2(.01f*c.SampleRate))));
		// if first element in bufferSizesSuggested is far from good size, insert goodSize as first element
		if (static_cast<float>(bufferSizesSuggested[0]) / goodSize < 0.666f || static_cast<float>(bufferSizesSuggested[0]) / goodSize > 1.5f)
		{
			bufferSizes.insert(bufferSizes.begin(), goodSize);
		}
		else // otherwise keep that element as first element
		{
			bufferSizes.insert(bufferSizes.begin() + 1, goodSize);
		}

		// return with first value that is a valid size
		for (auto& size : bufferSizes)
		{
			if (FFTReal::IsFFTSizeValid(size * CalculateNBuffersInFrame(c)))
			{
				c.BufferSize = size;
				return c;// return with first match
			}
		}
		return c; // failed to find valid BufferSize so use default
	}

	IIR2ndDF2Transposed FilterHighpass;
	FilterbankAnalysis Filterbank;
	FilterbankAnalysis FilterbankLoopback;
	EchoCancellerMomentum EchoCanceller;
	DetectVoiceActivation VAD;
	BeamformerAdaptive Beamformer;
	NoiseSuppression NoiseReduction;
	FilterbankSynthesis FilterbankInverse;

private:
	struct Coefficients
	{
		int NChannelsIn = 4; // Number of Input channels is NChannelsIn - NChannelsLoopback
		int NChannelsLoopback = 2;
		int BufferSize = 512;
		float SampleRate = 44100.f;
		float DelayDesired = .03f;
		AsynchronousBufferType AsynchronousBuffer = VARIABLE_SIZE;
	} C;

	struct Parameters
	{

	} P;

	struct Data
	{
		int LatencySamples;
		int NChannelsInput;
		int FFTSize, NBands, NBuffersInFrame;
		bool Activity;
		Eigen::ArrayXXf BufferTime;
		Eigen::ArrayXXcf BufferFreq, LoopbackFreq;
		void Reset() 
		{
			BufferTime.setZero();
			BufferFreq.setZero();
			LoopbackFreq.setZero();
			Activity = false;
		}
		bool InitializeMemory(const Coefficients& c)	
		{
			NChannelsInput = std::max(c.NChannelsIn - c.NChannelsLoopback, 0);
			NBuffersInFrame = CalculateNBuffersInFrame(c);
			LatencySamples = (NBuffersInFrame - 1) * c.BufferSize;
			FFTSize = NBuffersInFrame * c.BufferSize;
			NBands = FFTSize / 2 + 1;
			BufferTime.resize(c.BufferSize, NChannelsInput);
			BufferFreq.resize(NBands, NChannelsInput);
			LoopbackFreq.resize(NBands, NChannelsInput);
			return true;
		}
		size_t GetAllocatedMemorySize() const 
		{ 
			size_t size = BufferTime.GetAllocatedMemorySize();
			size += BufferFreq.GetAllocatedMemorySize();
			size += LoopbackFreq.GetAllocatedMemorySize();
			return size; 
		}
		void OnParameterChange(const Parameters& p, const Coefficients& c) { }
	} D;

	DEFINEMEMBERALGORITHMS(8, FilterHighpass, Filterbank, FilterbankLoopback, EchoCanceller, VAD, Beamformer, NoiseReduction, FilterbankInverse);

	auto InitializeMembers() 
	{
		const auto filterbankRate = C.SampleRate / C.BufferSize;
		// highpass filter
		auto cFH = FilterHighpass.GetCoefficients();
		cFH.NChannelsIn = D.NChannelsInput;
		cFH.BufferSize = C.BufferSize;
		cFH.SampleRate = C.SampleRate;
		auto flag = FilterHighpass.Initialize(cFH);
		// filterbank
		auto cFB = Filterbank.GetCoefficients();
		cFB.BufferSize = C.BufferSize;
		cFB.FrameSize = D.FFTSize;
		cFB.FFTSize = D.FFTSize;
		cFB.NChannels = D.NChannelsInput;
		flag &= Filterbank.Initialize(cFB);
		// filterbank loopback
		auto cFBL = FilterbankLoopback.GetCoefficients();
		cFBL.BufferSize = C.BufferSize;
		cFBL.FrameSize = D.FFTSize;
		cFBL.FFTSize = D.FFTSize;
		cFBL.NChannels = C.NChannelsLoopback;
		flag &= FilterbankLoopback.Initialize(cFBL);
		// echo canceller
		auto cEC = EchoCanceller.GetCoefficients();
		cEC.FilterbankRate = filterbankRate;
		cEC.FilterLengthMilliseconds = 100.f;
		cEC.NBands = D.NBands;
		cEC.nBuffersInAnalysisFrame = D.NBuffersInFrame;
		cEC.nBuffersInSynthesisFrame = D.NBuffersInFrame;
		cEC.NChannels = D.NChannelsInput;
		cEC.NChannelsLoopback = C.NChannelsLoopback;
		flag &= EchoCanceller.Initialize(cEC);
		// VAD
		auto cVAD = VAD.GetCoefficients();
		cVAD.FilterbankRate = filterbankRate;
		cVAD.NBands = D.NBands;
		cVAD.NChannels = D.NChannelsInput;
		flag &= VAD.Initialize(cVAD);
		// Beamformer
		auto cBF = Beamformer.GetCoefficients();
		cBF.FilterbankRate = filterbankRate;
		cBF.NBands = D.NBands;
		cBF.NChannels = D.NChannelsInput;
		flag &= Beamformer.Initialize(cBF);
		// NoiseReduction
		auto cNR = NoiseReduction.GetCoefficients();
		cNR.FilterbankRate = filterbankRate;
		cNR.NBands = D.NBands;
		cNR.NChannels = 1;
		flag &= NoiseReduction.Initialize(cNR);
		// Inverse filterbank
		auto cFBI = FilterbankInverse.GetCoefficients();
		cFBI.BufferSize = C.BufferSize;
		cFBI.FrameSize = D.FFTSize;
		cFBI.FFTSize = D.FFTSize;
		cFBI.NChannels = 1;
		flag &= FilterbankInverse.Initialize(cFBI);
		return flag;
	}

	void ProcessOn(Input xTime, Output yTime)
	{
		// high pass input (but not loopback)
		FilterHighpass.Process(xTime.leftCols(D.NChannelsInput), D.BufferTime);

		Filterbank.Process(D.BufferTime, D.BufferFreq);
		FilterbankLoopback.Process(xTime.middleCols(D.NChannelsInput, C.NChannelsLoopback), D.LoopbackFreq);

		EchoCanceller.Process({ D.BufferFreq, D.LoopbackFreq }, D.BufferFreq);

		VAD.Process(D.BufferFreq, D.Activity);

		Beamformer.Process({ D.BufferFreq, D.Activity }, D.BufferFreq.col(0));

		NoiseReduction.Process(D.BufferFreq.col(0), D.BufferFreq.col(0));

		FilterbankInverse.Process(D.BufferFreq.col(0), yTime);
	}

	void ProcessOff(Input xTime, Output yTime) { yTime = xTime.col(0); }
};