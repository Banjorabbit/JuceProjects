#pragma once
#include "../BaseClasses/PureCRTP.h"
#include "../FilterMinMax.h"
#include "../Utilities/ApproximationMath.h"

class DetectTonalPowerMin : public Base<DetectTonalPowerMin, I::Complex2D, O::Bool2D>
{
	friend Base<DetectTonalPowerMin, I::Complex2D, O::Bool2D>;

public:
	FilterMin FindMin;

private:
	struct Coefficients
	{
		int NBands = 2049;
		int NChannels = 2;
		float SampleRate = 44100.f;
		float FilterbankRate = 44100.f / 512.f;
		float WindowSizeFreqHz = 200.f;
	} C;

	struct Parameters
	{
		float TonalTConstant = 0.1f;
		float TonalThreshold = 2.f;
		float PowerCompression = 0.3f;
	} P;

	struct Data
	{
		Eigen::ArrayXXf PowerMinSmooth;
		float TonalLambda;
		void Reset()
		{
			PowerMinSmooth.setZero();
		}
		bool InitializeMemory(const Coefficients& c)
		{
			PowerMinSmooth.resize(c.NBands, c.NChannels);
			return true;
		}
		size_t GetAllocatedMemorySize() const
		{
			return PowerMinSmooth.GetAllocatedMemorySize();
		}
		void OnParameterChange(const Parameters& p, const Coefficients& c)
		{
			TonalLambda = 1 - expf(-1.f / std::max(1e-10f, c.FilterbankRate*p.TonalTConstant));
		}
	} D;

	DEFINEMEMBERALGORITHMS(1, FindMin);

	auto InitializeMembers()
	{
		auto cMin = FindMin.GetCoefficients();
		auto nFFT = (C.NBands - 1) * 2;
		auto length = C.WindowSizeFreqHz / C.SampleRate * nFFT;
		cMin.Length = static_cast<int>(std::round(length / 2.f)) * 2 + 1; // make odd
		cMin.NChannels = C.NChannels;
		bool flag = FindMin.Initialize(cMin);

		return flag;
	}

	void ProcessOn(Input xFreq, Output tonalDetection)
	{
		Eigen::ArrayXXf xPower = xFreq.abs2();
		for (auto channel = 0; channel < C.NChannels; channel++)
		{
			for (auto i = 0; i < C.NBands; i++)
			{
				xPower(i, channel) = fasterpow(xPower(i, channel), P.PowerCompression);
			}
		}
		Eigen::ArrayXXf xMin(C.NBands, C.NChannels);
		FindMin.Process(xPower, xMin);

		xPower -= xMin + D.PowerMinSmooth;
		D.PowerMinSmooth += D.TonalLambda * xPower;
		tonalDetection = D.PowerMinSmooth > xMin * P.TonalThreshold;
	}

	void ProcessOff(Input xFreq, Output tonalDetection) { tonalDetection.setConstant(false); }
};
