#pragma once
#include "../BaseClasses/PureCRTP.h"
#include "../Filterbank.h"
#include "IIR2ndDF.h"
#include "../FrequencyDomain/BeamformerAdaptive.h"
#include "../FrequencyDomain/Deconvolver.h"
#include "../FrequencyDomain/DetectVoiceActivation.h"
#include "../FrequencyDomain/GainCalculation.h"
#include "../FrequencyDomain/CriticalBands.h"

class SpatialNoiseReduction : public AsynchronousBase<SpatialNoiseReduction>
{
	friend Base<SpatialNoiseReduction>;

public:
	int GetNChannelsOut() const { return 1; }
	auto GetLatencySamples() const { return D.LatencySamples; }

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
			if (FFTReal::IsFFTSizeValid(size * 4))
			{
				c.BufferSize = size;
				return c;// return with first match
			}
		}
		return c; // failed to find valid BufferSize so use default
	}

	BeamformerAdaptive Beamformer;
	CriticalBands ConvertToCritBands;
	Deconvolver deconvolver;
	DetectVoiceActivation VAD;
	FilterbankAnalysis Filterbank;
	FilterbankSynthesis FilterbankInverse;
	GainCalculationSimple GainCalculation;
	IIR2ndDF2Transposed FilterHighpass;
	
private:
	struct Coefficients
	{
		int NChannelsIn = 3;
		int BufferSize = 128;
		float SampleRate = 16000.f;
		AsynchronousBufferType AsynchronousBuffer = VARIABLE_SIZE;
	} C;

	struct Parameters
	{
	} P;

	struct Data
	{
		int NBands, LatencySamples;
		void Reset() {	}
		bool InitializeMemory(const Coefficients& c) 
		{ 
			LatencySamples = (4 - 1) * c.BufferSize;
			NBands = c.BufferSize * 4 / 2 + 1; 
			return true; 
		}
		size_t GetAllocatedMemorySize() const { return 0; }
		void OnParameterChange(const Parameters& p, const Coefficients& c) { }
	} D;

	DEFINEMEMBERALGORITHMS(6, FilterHighpass, Filterbank, VAD, Beamformer, FilterbankInverse, deconvolver);

	auto InitializeMembers()
	{
		const auto filterbankRate = C.SampleRate / C.BufferSize;
		const auto fftSize = C.BufferSize * 4;
		// highpass filter
		auto sFH = FilterHighpass.GetSetup();
		sFH.Coefficients.NChannelsIn = C.NChannelsIn;
		sFH.Coefficients.BufferSize = C.BufferSize;
		sFH.Coefficients.SampleRate = C.SampleRate;
		sFH.Parameters.FilterType = sFH.Parameters.Highpass;
		sFH.Parameters.Frequency = 150;
		auto flag = FilterHighpass.Initialize(sFH);
		// filterbank
		auto cFB = Filterbank.GetCoefficients();
		cFB.BufferSize = C.BufferSize;
		cFB.FrameSize = fftSize;
		cFB.FFTSize = fftSize;
		cFB.NChannels = C.NChannelsIn;
		flag &= Filterbank.Initialize(cFB);
		// Deconvolver
		auto cDC = deconvolver.GetCoefficients();
		cDC.NBands = D.NBands;
		cDC.NChannels = C.NChannelsIn;
		cDC.NLow = 3 * (int)C.SampleRate / 16000;
		cDC.N = 20 * (int)C.SampleRate / 16000;
		flag &= deconvolver.Initialize(cDC);
		// VAD
		auto cVAD = VAD.GetCoefficients();
		cVAD.FilterbankRate = filterbankRate;
		cVAD.NBands = D.NBands;
		cVAD.NChannels = C.NChannelsIn;
		flag &= VAD.Initialize(cVAD);
		// Beamformer
		auto cBF = Beamformer.GetCoefficients();
		cBF.FilterbankRate = filterbankRate;
		cBF.NBands = D.NBands;
		cBF.NChannels = C.NChannelsIn;
		flag &= Beamformer.Initialize(cBF);
		auto pBF = Beamformer.GetParameters();
		pBF.EnableNoiseFilter = true;
		Beamformer.SetParameters(pBF);

		// critical bands (needs to be initialized before GainCalculation
		auto cCB = ConvertToCritBands.GetCoefficients();
		cCB.NBands = D.NBands;
		cCB.SampleRate = C.SampleRate;
		flag &= ConvertToCritBands.Initialize(cCB);

		auto sGC = GainCalculation.GetSetup();
		sGC.Coefficients.FilterbankRate = C.SampleRate / C.BufferSize;
		sGC.Coefficients.NBands = ConvertToCritBands.GetNBandsCritical();
		sGC.Coefficients.NChannels = 1;
		sGC.Parameters.MinimumdB = -10;
		flag &= GainCalculation.Initialize(sGC);

		// Inverse filterbank
		auto cFBI = FilterbankInverse.GetCoefficients();
		cFBI.BufferSize = C.BufferSize;
		cFBI.FrameSize = fftSize;
		cFBI.FFTSize = fftSize;
		cFBI.NChannels = 1;
		flag &= FilterbankInverse.Initialize(cFBI);

		return flag;
	}

	void ProcessOn(Input xTime, Output yTime)
	{
		Eigen::ArrayXXf bufferTime(C.BufferSize, C.NChannelsIn);
		FilterHighpass.Process(xTime, bufferTime);

		Eigen::ArrayXXcf bufferFreq(D.NBands, C.NChannelsIn);
		Filterbank.Process(bufferTime, bufferFreq);

		bool activity;
		VAD.Process(bufferFreq, activity);
		
		Eigen::ArrayXXcf deconvolved(D.NBands, C.NChannelsIn);
		deconvolver.Process({ bufferFreq, activity }, deconvolved);

		Beamformer.Process({ deconvolved, activity }, bufferFreq.col(0));
		Eigen::ArrayXf powerN = Beamformer.GetNoisePower();
		Eigen::ArrayXf powerX = bufferFreq.col(0).abs2();

		Eigen::ArrayXf critPower(ConvertToCritBands.GetNBandsCritical());
		ConvertToCritBands.Process(bufferFreq.col(0).abs2(), critPower);

		Eigen::ArrayXf critPowerNoise(ConvertToCritBands.GetNBandsCritical());
		ConvertToCritBands.Process(Beamformer.GetNoise().abs2(), critPowerNoise);

		Eigen::ArrayXf snr = ((critPower.max(1e-20f) / critPowerNoise.max(1e-20f).min(1e3f))).max(1.f);

		Eigen::ArrayXf critGain(ConvertToCritBands.GetNBandsCritical());
		GainCalculation.Process(snr, critGain);

		Eigen::ArrayXf gain(D.NBands);
		ConvertToCritBands.Inverse<float>(critGain, gain);
		bufferFreq.col(0) *= gain;

		FilterbankInverse.Process(bufferFreq.col(0), yTime);
	}

	void ProcessOff(Input xTime, Output yTime) { yTime.setZero(); }
};
