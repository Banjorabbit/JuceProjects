#pragma once
#include "../BaseClasses/PureCRTP.h"

// CriticalBands converts an FFT power spectrum into critical bands using the Bark scale.
//
// The class also has a public inverse function.
//
// author: Kristian Timm Andersen
class CriticalBands : public Base<CriticalBands>
{
	friend Base<CriticalBands>;

public:
	int GetNBandsCritical() const { return D.NBandsCritical; }
	Eigen::Map<const Eigen::ArrayXf> GetCenterFrequencies() const { return Eigen::Map<const Eigen::ArrayXf>(D.FrequencyCentersArrayHz, GetNBandsCritical()); }
	
	// template to allow both floats, complex and bools to be inverted
	template<typename T>
	void Inverse(I::InArray2D<T> xPower, O::OutArray2D<T> yPower)
	{
		for (auto channel = 0; channel < xPower.cols(); channel++)
		{
			yPower.block(0, channel, D.IndexStart(0), 1).setConstant(xPower(0, channel));
			for (auto i = 0; i < D.NBandsCritical; i++)
			{
				yPower.block(D.IndexStart(i), channel, D.NSumBands(i), 1).setConstant(xPower(i, channel));
			}
			int indexEnd = D.IndexStart(D.NBandsCritical - 1) + D.NSumBands(D.NBandsCritical - 1);
			yPower.block(indexEnd, channel, C.NBands - indexEnd, 1).setConstant(xPower(D.NBandsCritical - 1, channel));
		}
	}

private:
	struct Coefficients
	{
		float SampleRate = 44.1e3f;
		int NBands = 257;
	} C;

	struct Parameters { } P;

	struct Data
	{
		const static float FrequencyCornersArrayHz[25]; // defined in .cpp file
		const static float FrequencyCentersArrayHz[24]; // defined in .cpp file
		int NBandsCritical = 24; // defined here so GetNBandsCritical() works before initialization
		Eigen::ArrayXi IndexStart, NSumBands;
		void Reset() { }
		bool InitializeMemory(const Coefficients& c)
		{
			// calculate number of critical bands
			Eigen::Map<const Eigen::ArrayXf> FrequencyCornersHz(FrequencyCornersArrayHz, 25);
			NBandsCritical = (FrequencyCornersHz < c.SampleRate / 2).cast<int>().sum() - 1;
			// convert Frequency corners to FFT indices
			auto NFFT = (c.NBands - 1) * 2; // number of FFT points
			Eigen::ArrayXi FrequencyCornersIndex = (FrequencyCornersHz * NFFT / c.SampleRate).round().cast<int>();
			IndexStart.resize(NBandsCritical);
			IndexStart = FrequencyCornersIndex.head(NBandsCritical);
			NSumBands.resize(NBandsCritical);
			NSumBands = (FrequencyCornersIndex.segment(1, NBandsCritical) - IndexStart).cwiseMax(1);
			return true;
		}
		size_t GetAllocatedMemorySize() const { return IndexStart.GetAllocatedMemorySize() + NSumBands.GetAllocatedMemorySize(); }
		void OnParameterChange(const Parameters& p, const Coefficients& c) { }
	} D;

	void ProcessOn(Input xPower, Output yPower)
	{
		for (auto channel = 0; channel < xPower.cols(); channel++)
		{
			for (auto i = 0; i < D.NBandsCritical; i++)
			{
				// mean in a critical band
				yPower(i, channel) = xPower.block(D.IndexStart(i), channel, D.NSumBands(i), 1).mean();
			}
		}

	}

	void ProcessOff(Input xPower, Output yPower) { yPower.setZero(); }
};

