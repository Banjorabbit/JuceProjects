#pragma once
#include "../BaseClasses/PureCRTP.h"
#include "../FrequencyDomain/InterpolationTemporal.h"
#include "../Filterbank.h"
#include "../InterpolationCubic.h"

// TODO: Support for more than 1 channel
// TODO: Investigate time-varying P.StretchFactor. Consider replacing D.OutBuffer with circularBuffer (CircBuffer.h)
class PitchShift : public AsynchronousBase<PitchShift>
{
	friend Base<PitchShift>;

public:
	InterpolationTemporal InterpolateSpectrogram;
	FilterbankAnalysis Filterbank;
	FilterbankSynthesis FilterbankInverse;
	InterpolationCubic InterpolateTime;

	int GetLatencySamples() const { return D.LatencySamples; }

	auto GetAsynchronousCoefficients(const std::vector<int> bufferSizesSuggested) const
	{
		Coefficients c = C;
		auto bufferSizes = bufferSizesSuggested;
			   
		// start suggestion should be size 2^x close to 3ms
		int goodSize = static_cast<int>(powf(2, std::round(std::log2(.003f*c.SampleRate))));
		if (static_cast<float>(bufferSizesSuggested[0]) / goodSize < 0.666f || static_cast<float>(bufferSizesSuggested[0]) / goodSize > 1.5f)
		{
			bufferSizes.insert(bufferSizes.begin(), goodSize);
		}
		else
		{
			bufferSizes.insert(bufferSizes.begin() + 1, goodSize);
		}

		// return with first value that is a valid size
		for (auto& size : bufferSizes)
		{
			if (FFTReal::IsFFTSizeValid(size * 4))
			{
				c.BufferSize = size;
				return c; // return with first match
			}
		}
		return c; // failed to find valid BufferSize so use default
	}

private:
	struct Coefficients
	{
		int BufferSize = 128;
		ConstrainedType<int> NChannelsIn = { 1, 1, 1 };
		float SampleRate = 44.1e3f;
		AsynchronousBufferType AsynchronousBuffer = VARIABLE_SIZE;
	} C;

	struct Parameters
	{
		float StretchFactor = 1.f;
	} P;

	struct Data
	{
		int LatencySamples;
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
			IndexOut = static_cast<int>(OutBuffer.size()) / 2 + 1;
			LatencySamples = IndexOut;
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
		sF.Coefficients.NChannels = C.NChannelsIn;
		sF.Coefficients.FrameSize = D.FFTSize;
		sF.Coefficients.FFTSize = D.FFTSize;
		sF.Coefficients.BufferSize = C.BufferSize;
		sF.Parameters.WindowType = sF.Parameters.HannWindow;
		flag &= Filterbank.Initialize(sF);

		auto sFI = FilterbankInverse.GetSetup();
		sFI.Coefficients.NChannels = C.NChannelsIn;
		sFI.Coefficients.FrameSize = D.FFTSize;
		sFI.Coefficients.FFTSize = D.FFTSize;
		sFI.Coefficients.BufferSize = C.BufferSize;
		sFI.Parameters.WindowType = sFI.Parameters.HannWindow;
		flag &= FilterbankInverse.Initialize(sFI);

		flag &= InterpolateTime.Initialize();
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
		D.OutBuffer.head(3 * C.BufferSize) = D.OutBuffer.tail(3 * C.BufferSize);
		D.OutBuffer.tail(C.BufferSize).setZero();
		D.IndexOut -= C.BufferSize;
	}

	void ProcessOff(Input xTime, Output yTime) { yTime = xTime; }
};
