#pragma once
#include "../BaseClasses/PureCRTP.h"
#include "PitchShift.h"
#include "../Utilities//ApproximationMath.h"
#include "../Filterbank.h"
#include "../FrequencyDomain/DetectTransient.h"

class PitchShiftAdaptiveResolution : public Base<PitchShiftAdaptiveResolution>
{
	friend Base<PitchShiftAdaptiveResolution>;

public:
	PitchShift PitchShiftShortFrame;
	PitchShift PitchShiftLongFrame;
	VectorAlgo<FilterbankAnalysis> Filterbank;
	FilterbankSynthesis FilterbankInverse;
	DetectTransient TransientDetection;

	// Latency is Long Pitchshift delay + 1 filterbank delay
	int GetLatencySamples() const { return PitchShiftLongFrame.GetLatencySamples() + C.BufferSize * 3; }
private:

	struct Coefficients 
	{
		float SampleRate = 44.1e3f;
		int BufferSize = 128;
		int LongFrameFactor = 5;
	} C;
	
	struct Parameters 
	{
		float StretchFactor = 1.f;
	} P;

	struct Data 
	{
		Eigen::ArrayXf xTimeLong, yTimeLong, yTimeShort, PhaseOld;
		decltype(PitchShiftShortFrame.GetParameters()) pShort;
		decltype(PitchShiftLongFrame.GetParameters()) pLong;
		int IndexLong, FFTSize, NBands;
		void Reset() 
		{
			xTimeLong.setZero();
			yTimeLong.setZero();
			yTimeShort.setZero();
			PhaseOld.setZero();
			IndexLong = 0;
		}
		bool InitializeMemory(const Coefficients& c)
		{
			xTimeLong.resize(c.BufferSize * c.LongFrameFactor);
			yTimeLong.resize(c.BufferSize * c.LongFrameFactor);
			yTimeShort.resize(c.BufferSize * c.LongFrameFactor);
			FFTSize = 4 * c.BufferSize;
			NBands = FFTSize / 2 + 1;
			PhaseOld.resize(NBands);
			return true;
		}
		size_t GetAllocatedMemorySize() const 
		{
			size_t size = xTimeLong.GetAllocatedMemorySize();
			size += yTimeLong.GetAllocatedMemorySize();
			size += yTimeShort.GetAllocatedMemorySize();
			size += PhaseOld.GetAllocatedMemorySize();
			return size;
		}
		void OnParameterChange(const Parameters& p, const Coefficients& c) 
		{
			pShort.StretchFactor = p.StretchFactor;
			pLong.StretchFactor = p.StretchFactor;
		}
	} D;

	DEFINEMEMBERALGORITHMS(5, PitchShiftShortFrame, PitchShiftLongFrame, Filterbank, FilterbankInverse, TransientDetection);

	auto InitializeMembers()
	{
		auto sShort = PitchShiftShortFrame.GetSetup();
		sShort.Coefficients.BufferSize = C.BufferSize;
		bool flag = PitchShiftShortFrame.Initialize(sShort);

		auto sLong = PitchShiftLongFrame.GetSetup();
		sLong.Coefficients.BufferSize = C.BufferSize * C.LongFrameFactor;
		flag &= PitchShiftLongFrame.Initialize(sLong);

		Filterbank.resize(2);
		auto sFilterbank = Filterbank[0].GetSetup();
		sFilterbank.Coefficients.BufferSize = C.BufferSize;
		sFilterbank.Coefficients.FFTSize = D.FFTSize;
		sFilterbank.Coefficients.FrameSize = D.FFTSize;
		sFilterbank.Coefficients.NChannels = 1; // change when more channels are supported
		sFilterbank.Parameters.WindowType = sFilterbank.Parameters.HannWindow;
		flag &= Filterbank.Initialize(sFilterbank);

		auto sFilterbankInverse = FilterbankInverse.GetSetup();
		sFilterbankInverse.Coefficients.BufferSize = C.BufferSize;
		sFilterbankInverse.Coefficients.FFTSize = D.FFTSize;
		sFilterbankInverse.Coefficients.FrameSize = D.FFTSize;
		sFilterbankInverse.Coefficients.NChannels = 1; // change when more channels are supported
		sFilterbankInverse.Parameters.WindowType = sFilterbankInverse.Parameters.HannWindow;
		flag &= FilterbankInverse.Initialize(sFilterbankInverse);

		auto cTD = TransientDetection.GetCoefficients();
		cTD.FilterbankRate = C.SampleRate / C.BufferSize;
		cTD.NBands = D.NBands;
		cTD.NChannels = 1; // change when more channels are supported
		cTD.SampleRate = C.SampleRate;
		flag &= TransientDetection.Initialize(cTD);

		return flag;
	}

	void ProcessOn(Input xTime, Output yTime) 
	{
		// process Long and Short PitchShift
		PitchShiftShortFrame.Process(xTime, D.yTimeShort.segment(D.IndexLong * C.BufferSize, C.BufferSize));
		D.xTimeLong.segment(D.IndexLong * C.BufferSize, C.BufferSize) = xTime;
		D.IndexLong++;
		if (D.IndexLong == C.LongFrameFactor) 
		{ 
			D.IndexLong = 0;
			PitchShiftLongFrame.Process(D.xTimeLong, D.yTimeLong); 
			// set after processing to ensure Long and Short are syncronized
			PitchShiftShortFrame.SetParameters(D.pShort);
			PitchShiftLongFrame.SetParameters(D.pLong);
		}

		// Detect transience and determine which Pitchshift to use in each time/frequency bin
		Eigen::ArrayXXcf yFreq(D.NBands, 2);
		Filterbank[0].Process(D.yTimeShort.segment(D.IndexLong * C.BufferSize, C.BufferSize), yFreq.col(0));
		Filterbank[1].Process(D.yTimeLong.segment(D.IndexLong * C.BufferSize, C.BufferSize), yFreq.col(1));
		
		Eigen::Array<bool, Eigen::Dynamic, 1> transientFlag(D.NBands), transientFlagCritical(TransientDetection.ConvertBands.GetNBandsCritical());
		TransientDetection.Process(yFreq.col(0).abs2(), transientFlagCritical);
		TransientDetection.ConvertBands.Inverse<bool>(transientFlagCritical, transientFlag);
		
		for (auto i = 0; i < D.NBands; i++) 
		{ 
			if (transientFlag(i)) 
			{ 
				float phase = std::arg(yFreq(i, 0)) + D.PhaseOld(i);
				float mag = std::abs(yFreq(i, 0));
				yFreq.real()(i, 1) = mag * std::cos(phase);
				yFreq.imag()(i, 1) = mag * std::sin(phase);
			} 
		}
		D.PhaseOld = yFreq.col(1).arg() - yFreq.col(0).arg();
		for (auto i = 0; i < D.NBands; i++)
		{
			// wrap to +- PI (0 Hz always has Phase=0.f)
			if (D.PhaseOld(i) > static_cast<const float>(PI)) { D.PhaseOld(i) -= 2.f*static_cast<const float>(PI); }
			else if (D.PhaseOld(i) < -static_cast<const float>(PI)) { D.PhaseOld(i) += 2.f*static_cast<const float>(PI); }
		}

		FilterbankInverse.Process(yFreq.col(1), yTime);
	}

	void ProcessOff(Input xTime, Output yTime) { yTime.setZero(); }
};


class PitchShiftAdaptiveResolutionStreaming : public AsynchronousStreaming<PitchShiftAdaptiveResolutionStreaming, PitchShiftAdaptiveResolution>
{
	friend AsynchronousStreaming<PitchShiftAdaptiveResolutionStreaming, PitchShiftAdaptiveResolution>;

	// latency is equal to reset-value of D.IndexOut 
	int CalculateLatencySamples() const { return Algo.GetLatencySamples(); }
	int CalculateNChannelsOut(const int nChannels) const { return nChannels; }

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