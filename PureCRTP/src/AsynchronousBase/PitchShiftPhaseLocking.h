#pragma once
#include "../BaseClasses/PureCRTP.h"
#include "../FrequencyDomain/InterpolationTemporal.h"
#include "../FrequencyDomain/DetectTonal.h"
#include "../Filterbank.h"
#include "../InterpolationCubic.h"

// TODO: Support for more than 1 channel
// TODO: Investigate time-varying P.StretchFactor. Consider replacing D.OutBuffer with circularBuffer (CircBuffer.h)
class PitchShiftPhaseLocking : public AsynchronousBase<PitchShiftPhaseLocking>
{
	friend Base<PitchShiftPhaseLocking>;

public:
	DetectTonal TonalDetector;
	InterpolationTemporal InterpolateSpectrogram;
	FilterbankAnalysis Filterbank;
	FilterbankSynthesis FilterbankInverse;
	FilterbankSynthesis FilterbankInverseResidual;
	InterpolationCubic InterpolateTime;

	int GetLatencySamples() const { return D.LatencySamples; }

	int GetNChannelsOut() const
	{
		return C.NChannelsIn * 2;
	}

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
		int BufferSize = 512;
		int FFTSize = 4096;
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
			IndexOut = static_cast<int>(OutBuffer.size()) / 2 + 1;
			Energy.setZero();
			Phase.setZero();
			FractionalDelay = 2.f; // make first delay = 1
			InverseDelay = 2.f;
			OutBuffer.setZero();
			SynthesisBuffer.setZero();
		}
		bool InitializeMemory(const Coefficients& c)
		{
			NBands = c.FFTSize / 2 + 1;
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

	DEFINEMEMBERALGORITHMS(6, TonalDetector, InterpolateSpectrogram, Filterbank, FilterbankInverse, FilterbankInverseResidual, InterpolateTime)

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

		auto cI = InterpolateSpectrogram.GetCoefficients();
		cI.NBands = D.NBands;
		flag &= InterpolateSpectrogram.Initialize(cI);

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

		auto sFIR = FilterbankInverseResidual.GetSetup();
		sFIR.Coefficients.NChannels = C.NChannelsIn;
		sFIR.Coefficients.FrameSize = C.FFTSize;
		sFIR.Coefficients.FFTSize = C.FFTSize;
		sFIR.Coefficients.BufferSize = C.BufferSize;
		sFIR.Parameters.WindowType = sFIR.Parameters.HannWindow;
		sFIR.Parameters.Gain = static_cast<float>(4 * C.BufferSize) / C.FFTSize;
		flag &= FilterbankInverseResidual.Initialize(sFIR);

		flag &= InterpolateTime.Initialize();
		return flag;
	}

	void ProcessOn(Input xTime, Output yTime)
	{
		Eigen::ArrayXcf xFreq(D.NBands);
		Filterbank.Process(xTime, xFreq);

		Eigen::Array<bool, Eigen::Dynamic, 1> decisionTonal(D.NBands);
		TonalDetector.Process(xFreq, decisionTonal);

		// find tonal part
		Eigen::ArrayXcf xTonal = xFreq * decisionTonal.cast<float>();
		// synthesize residual (not tonal)
		FilterbankInverseResidual.Process(xFreq * (1.f - decisionTonal.cast<float>()), yTime.col(1));

		D.Energy.col(0) = D.Energy.col(1);
		D.Energy.col(1) = xTonal.abs2();
		D.Phase.col(0) = D.Phase.col(1);
		D.Phase.col(1) = xTonal.arg();
		D.FractionalDelay -= 1.f; // Decrease by 1 since we have 1 new frame
		for (D.FractionalDelay; D.FractionalDelay < 1.f; D.FractionalDelay += P.StretchFactor)
		{
			InterpolateSpectrogram.Process({ D.Energy, D.Phase, D.FractionalDelay }, xTonal);
			D.SynthesisBuffer.head(3) = D.SynthesisBuffer.tail(3);
			FilterbankInverse.Process(xTonal, D.SynthesisBuffer.tail(C.BufferSize));
			for (D.InverseDelay; D.InverseDelay < C.BufferSize; D.InverseDelay += D.InverseFactor)
			{
				auto iDelay = static_cast<int>(D.InverseDelay);
				Eigen::Array<float, 1, 1> fractionalDelay;
				fractionalDelay(0) = D.InverseDelay - static_cast<float>(iDelay);
				InterpolateTime.Process({ Eigen::Map<Eigen::Array4f>(&D.SynthesisBuffer(iDelay)), fractionalDelay }, D.OutBuffer.segment(D.IndexOut, 1));
				D.IndexOut++;
			}
			D.InverseDelay -= C.BufferSize;
		}
		yTime.col(0) = D.OutBuffer.head(C.BufferSize);
		D.OutBuffer.head(3 * C.BufferSize) = D.OutBuffer.tail(3 * C.BufferSize);
		D.OutBuffer.tail(C.BufferSize).setZero();
		D.IndexOut -= C.BufferSize;

		
	}

	void ProcessOff(Input xTime, Output yTime) { yTime = xTime.replicate(1,2); }
};
