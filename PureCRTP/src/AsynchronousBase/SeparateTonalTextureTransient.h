#pragma once
#include "../BaseClasses/PureCRTP.h"
#include "SeparateTransient.h"
#include "SeparateTonal.h"
#include "../CircBuffer.h"

class SeparateTonalTextureTransient : public AsynchronousBase<SeparateTonalTextureTransient>
{
	friend Base<SeparateTonalTextureTransient>;

public:
	SeparateTonal TonalSeparator;
	SeparateTransient TransientSeparator;
	CircBuffer DelayLine;

	int GetLatencySamples() const { return TonalSeparator.GetLatencySamples() + TransientSeparator.GetLatencySamples(); }

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

	struct Parameters 
	{
		enum SelectOutput { TONAL, TEXTURE, TRANSIENT, SUM};
		SelectOutput OutputSelector = TONAL;
	} P;

	struct Data
	{
		int NBands;
		void Reset() { }
		bool InitializeMemory(const Coefficients& c)
		{
			NBands = c.FFTSize / 2 + 1;
			return true;
		}
		size_t GetAllocatedMemorySize() const { return 0; }
		void OnParameterChange(const Parameters& p, const Coefficients& c) { }
	} D;

	DEFINEMEMBERALGORITHMS(3, TonalSeparator, TransientSeparator, DelayLine)

	auto InitializeMembers()
	{
		auto sTonal = TonalSeparator.GetSetup();
		sTonal.Coefficients.BufferSize = C.BufferSize;
		sTonal.Coefficients.FrameSize = C.FFTSize;
		sTonal.Coefficients.NChannelsIn = C.NChannelsIn;
		sTonal.Coefficients.SampleRate = C.SampleRate;
		auto flag = TonalSeparator.Initialize(sTonal);

		auto sTransient = TransientSeparator.GetSetup();
		sTransient.Coefficients.BufferSize = C.BufferSize;
		sTransient.Coefficients.NChannelsIn = C.NChannelsIn;
		sTransient.Coefficients.SampleRate = C.SampleRate;
		flag &= TransientSeparator.Initialize(sTransient);

		auto cDelayLine = DelayLine.GetCoefficients();
		cDelayLine.DelayLength = TonalSeparator.GetLatencySamples(); // TonalSeparator must be initialized before this call
		cDelayLine.NChannels = C.NChannelsIn;
		flag &= DelayLine.Initialize(cDelayLine);

		return flag;
	}

	void ProcessOn(Input xTime, Output yTime)
	{
		Eigen::ArrayXXf xTonal(C.BufferSize, C.NChannelsIn), xTexture(C.BufferSize, C.NChannelsIn), xTransient(C.BufferSize, C.NChannelsIn);
		TonalSeparator.Process(xTime, xTonal);
		DelayLine.Process(xTime, xTexture);
		xTexture -= xTonal;
		TransientSeparator.Process(xTexture, xTransient);
		//xTexture -= xTransient;		
		switch (P.OutputSelector)
		{
		case Parameters::TONAL: yTime = xTonal; break;
		case Parameters::TEXTURE: yTime = xTexture; break;
		case Parameters::TRANSIENT: yTime = xTransient; break;
		case Parameters::SUM: // go to default
		default: yTime = xTonal + xTexture + xTransient;
		}
	}

	void ProcessOff(Input xTime, Output yTime) { yTime = xTime; }
};