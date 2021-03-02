#pragma once
#include "../BaseClasses/PureCRTP.h"
#include "../FFT.h"
#include "../Utilities/ApproximationMath.h"

// Detect transients.
//
// author: Kristian Timm Andersen
class DetectTransientSmooth : public Base<DetectTransientSmooth, I::Complex2D, O::Bool2D>
{
	friend Base<DetectTransientSmooth, I::Complex2D, O::Bool2D>;

public:
	FFTReal FFT;

private:
	struct Coefficients
	{
		int NBands = 257;
		float SampleRate = 44.1e3f;
		float FilterbankRate = 44.1e3f / 256;
		int NChannels = 2;
		float WindowSizeFreqHz = 800.f;
	} C;

	struct Parameters
	{
		float TConstant = 0.1f;
		float Threshold = 1.5f;
		float PowerCompression = 0.3f;
	} P;

	struct Data
	{
		float Lambda;
		Eigen::ArrayXcf HFilter;
		Eigen::ArrayXXf PowerSmooth;
		void Reset()
		{
			PowerSmooth.setZero();
		}
		bool InitializeMemory(const Coefficients& c)
		{
			HFilter.resize(c.NBands);
			PowerSmooth.resize(c.NBands, c.NChannels);
			return true;
		}
		size_t GetAllocatedMemorySize() const
		{
			size_t size = HFilter.GetAllocatedMemorySize();
			size += PowerSmooth.GetAllocatedMemorySize();
			return size;
		}
		void OnParameterChange(const Parameters& p, const Coefficients& c)
		{
			Lambda = 1.f - expf(-1.f / std::max(1e-10f, c.FilterbankRate * p.TConstant));
		}
	} D;

	DEFINEMEMBERALGORITHMS(1, FFT);

	auto InitializeMembers()
	{
		auto nFFT = (C.NBands - 1) * 2;
		auto length = C.WindowSizeFreqHz / C.SampleRate * nFFT;
		
		auto cFFT = FFT.GetCoefficients();
		cFFT.FFTSize = nFFT;
		auto flag = FFT.Initialize(cFFT);

		// initialize D.HFilter using FFT
		Eigen::ArrayXf h = Eigen::ArrayXf::Zero(nFFT);
		h.head(length / 2) = 2.f / length; // averaging filter
		FFT.Process(h, D.HFilter);
		D.HFilter = D.HFilter.abs2(); // square averaging filter in freq domain

		return flag;
	}

	void ProcessOn(Input xFreq, Output activityFlag)
	{
		for (auto channel = 0; channel < C.NChannels; channel++)
		{
			Eigen::ArrayXf xPower = xFreq.col(channel).abs2();
			for (auto i = 0; i < C.NBands; i++)
			{
				xPower(i) = fasterpow(xPower(i), P.PowerCompression);
			}
			Eigen::ArrayXf xPowerMirror(2 * (C.NBands - 1));
			xPowerMirror.head(C.NBands) = xPower;
			xPowerMirror.tail(C.NBands - 2) = xPower.segment(1, C.NBands - 2).colwise().reverse();
			Eigen::ArrayXcf convFilter(C.NBands);
			FFT.Process(xPowerMirror, convFilter);
			convFilter *= D.HFilter;
			FFT.Inverse(convFilter, xPowerMirror);
			D.PowerSmooth.col(channel) += D.Lambda * (xPowerMirror.head(C.NBands) - D.PowerSmooth.col(channel));
			activityFlag.col(channel) = xPowerMirror.head(C.NBands) > D.PowerSmooth.col(channel) * P.Threshold;
		}
	}

	void ProcessOff(Input xFreq, Output activityFlag) { activityFlag.setConstant(false); }
};