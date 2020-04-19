#pragma once
#include "../BaseClasses/PureCRTP.h"
#include "../FFT.h"

// Calculate minimum phase spectrum from magnitude spectrum
//
// author: Kristian Timm Andersen

class MinPhaseSpectrum : public BaseFrequencyDomain<MinPhaseSpectrum, I::Real2D>
{
	friend BaseFrequencyDomain<MinPhaseSpectrum, I::Real2D>;

public:
	FFTReal FFT;

private:
	struct Coefficients 
	{
		int NBands = 257;
	} C;

	struct Parameters { } P;

	struct Data 
	{
		void Reset() {}
		bool InitializeMemory(const Coefficients& c) { return true; }
		size_t GetAllocatedMemorySize() const {	return 0; }
		void OnParameterChange(const Parameters& p, const Coefficients& c) {}
	} D;

	DEFINEMEMBERALGORITHMS(1, FFT);

	auto InitializeMembers()
	{
		auto setup = FFT.GetSetup();
		setup.Coefficients.FFTSize = (C.NBands - 1) * 2;
		return FFT.Initialize(setup);
	}

	void ProcessOn(Input xMag, Output yFreq) 
	{
		using namespace std::complex_literals;

		for (auto channel = 0; channel < xMag.cols(); channel++)
		{
			// calculate cepstrum
			Eigen::ArrayXcf xLog = xMag.col(channel).max(1e-5f).log().cast<std::complex<float>>();
			Eigen::ArrayXf xCepstrum((C.NBands - 1) * 2);
			FFT.Inverse(xLog, xCepstrum);
			// fold
			xCepstrum.segment(1, C.NBands - 2) += xCepstrum.segment(C.NBands, C.NBands - 2).colwise().reverse();
			xCepstrum.segment(C.NBands, C.NBands - 2) = 0.f;
			// convert back
			FFT.Process(xCepstrum, xLog);
			yFreq.col(channel) = xLog.real().exp() * (xLog.imag().cos() + 1.if*xLog.imag().sin()); //complex exponential exp(x+iy) = exp(x)*(cos(y)+i*sin(y))
		}
	}

	void ProcessOff(Input xMag, Output yFreq) { yFreq = xMag.cast<std::complex<float>>(); }
};