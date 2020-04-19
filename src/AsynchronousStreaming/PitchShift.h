#pragma once
#include "../BaseClasses/PureCRTP.h"
#include "../FrequencyDomain/InterpolationTemporal.h"
#include "../Filterbank.h"
#include "../InterpolationCubic.h"

// TODO: Support for more than 1 channel
// TODO: Investigate time-varying P.StretchFactor. Will D.Outbuffer overflow? Consider replacing D.OutBuffer with circularBuffer (CircBuffer.h)
class PitchShift : public Base<PitchShift>
{
	friend Base<PitchShift>;

public:
	InterpolationTemporal InterpolateSpectrogram;
	FilterbankAnalysis Filterbank;
	FilterbankSynthesis FilterbankInverse;
	InterpolationCubic InterpolateTime;

private:
	struct Coefficients
	{
		int BufferSize = 128;
		ConstrainedType<int> NChannels = { 1, 1, 1 };
	} C;

	struct Parameters
	{
		float StretchFactor = .8f;
	} P;

	struct Data
	{
		int FFTSize, NBands, IndexOut;
		float FractionalDelay, InverseFactor, InverseDelay;
		Eigen::ArrayX2f Energy, Phase;
		Eigen::ArrayXf OutBuffer, SynthesisBuffer;
		void Reset()
		{
			IndexOut = static_cast<int>(OutBuffer.size())/2 + 1;
			Energy.setZero();
			Phase.setZero();
			FractionalDelay = 2.f; // make first delay = 1
			InverseDelay = 2.f;
			OutBuffer.setZero();
			SynthesisBuffer.setZero();
		}
		bool InitializeMemory(const Coefficients& c)
		{
			FFTSize = c.BufferSize * 4;
			NBands = FFTSize / 2 + 1;
			Energy.resize(NBands, 2);
			Phase.resize(NBands, 2);
			OutBuffer.resize(c.BufferSize * 4); // make big to accommodate changing P.stretchfactors
			SynthesisBuffer.resize(c.BufferSize + 3);
			return true;
		}
		size_t GetAllocatedMemorySize() const
		{
			size_t size = Energy.GetAllocatedMemorySize();
			size += Phase.GetAllocatedMemorySize();
			size += OutBuffer.GetAllocatedMemorySize();
			size += SynthesisBuffer.GetAllocatedMemorySize();
			return size;
		}
		void OnParameterChange(const Parameters& p, const Coefficients& c) 
		{
			InverseFactor = 1.f / p.StretchFactor;
		}
	} D;

	DEFINEMEMBERALGORITHMS(4, InterpolateSpectrogram, Filterbank, FilterbankInverse, InterpolateTime)

	auto InitializeMembers()
	{
		auto cI = InterpolateSpectrogram.GetCoefficients();
		cI.NBands = D.NBands;
		auto flag = InterpolateSpectrogram.Initialize(cI);

		auto sF = Filterbank.GetSetup();
		sF.Coefficients.NChannels = 1;
		sF.Coefficients.FrameSize = D.FFTSize;
		sF.Coefficients.FFTSize = D.FFTSize;
		sF.Coefficients.BufferSize = C.BufferSize;
		sF.Parameters.WindowType = sF.Parameters.HannWindow;
		flag &= Filterbank.Initialize(sF);

		auto sFI = FilterbankInverse.GetSetup();
		sFI.Coefficients.NChannels = 1;
		sFI.Coefficients.FrameSize = D.FFTSize;
		sFI.Coefficients.FFTSize = D.FFTSize;
		sFI.Coefficients.BufferSize = C.BufferSize;
		sFI.Parameters.WindowType = sFI.Parameters.HannWindow;
		flag &= FilterbankInverse.Initialize(sFI);
		return flag;
	}

	void ProcessOn(Input xTime, Output yTime)
	{
		Eigen::ArrayXcf xFreq(D.NBands);
		Filterbank.Process(xTime, xFreq);
		D.Energy.col(0) = D.Energy.col(1);
		D.Energy.col(1) = xFreq.abs2();
		D.Phase.col(0) = D.Phase.col(1);
		D.Phase.col(1) = xFreq.arg();
		D.FractionalDelay -= 1.f; // Decrease by 1 since we have 1 new frame
		for (D.FractionalDelay; D.FractionalDelay < 1.f; D.FractionalDelay += P.StretchFactor)
		{
			InterpolateSpectrogram.Process({ D.Energy, D.Phase, D.FractionalDelay }, xFreq);
			D.SynthesisBuffer.head(3) = D.SynthesisBuffer.tail(3);
			FilterbankInverse.Process(xFreq, D.SynthesisBuffer.tail(C.BufferSize));
			for (D.InverseDelay; D.InverseDelay < C.BufferSize; D.InverseDelay += D.InverseFactor)
			{
				auto iDelay = static_cast<int>(D.InverseDelay);
				Eigen::Array<float, 1, 1> fractionalDelay;
				fractionalDelay(0) =  D.InverseDelay - static_cast<float>(iDelay);
				InterpolateTime.Process({ Eigen::Map<Eigen::Array4f>(&D.SynthesisBuffer(iDelay)), fractionalDelay }, D.OutBuffer.segment(D.IndexOut,1));
				D.IndexOut++;
			}
			D.InverseDelay -= C.BufferSize;
		}
		yTime = D.OutBuffer.head(C.BufferSize);
		D.OutBuffer.head(C.BufferSize) = D.OutBuffer.tail(C.BufferSize);
		D.IndexOut -= C.BufferSize;
	}

	void ProcessOff(Input xTime, Output yTime) { yTime = xTime; }
};

class PitchShiftStreaming : public AsynchronousStreaming<PitchShiftStreaming, PitchShift>
{
	friend AsynchronousStreaming<PitchShiftStreaming, PitchShift>;

	// latency is equal to reset-value of D.IndexOut 
	int GetLatencySamples(const decltype(Algo.GetCoefficients())& c) const { return c.BufferSize * 2 + 1; }
	int GetNChannelsOut(const int nChannels) const { return nChannels; }

	int UpdateCoefficients(decltype(Algo.GetCoefficients())& c, const float sampleRate, const int nChannels, const int bufferSizeExpected, std::vector<int> bufferSizesSuggested) const
	{
		if (nChannels > 1) { return -1; }

		// start suggestion should be size 2^x close to 3ms
		int goodSize = static_cast<int>(powf(2, std::round(std::log2(.003f*sampleRate))));
		if (static_cast<float>(bufferSizesSuggested[0]) / goodSize < 0.666f || static_cast<float>(bufferSizesSuggested[0]) / goodSize > 1.5f)
		{
			bufferSizesSuggested.insert(bufferSizesSuggested.begin(), goodSize);
		}
		else
		{
			bufferSizesSuggested.insert(bufferSizesSuggested.begin() + 1, goodSize);
		}
		bufferSizesSuggested.push_back(c.BufferSize);

		// return with first value that is a valid size
		for (auto& size : bufferSizesSuggested)
		{
			if (FFTReal::IsFFTSizeValid(size * 4))
			{
				c.BufferSize = size;
				return size;
			}
		}
		// failed to set up correctly
		return -1;
	}
};