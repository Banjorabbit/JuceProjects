#pragma once
#include "../BaseClasses/PureCRTP.h"
#include "../FilterMinMax.h"
#include "../FFT.h"
#include "../Utilities/ApproximationMath.h"

class DetectTonal : public Base<DetectTonal, I::Complex2D, O::Bool2D>
{
	friend Base<DetectTonal, I::Complex2D, O::Bool2D>;

public:
	FilterMin FindMin;
	FFTReal FFT;

private:
	struct Coefficients
	{
		int NBands = 1025;
		int NChannels = 2;
		float SampleRate = 44100.f;
		float FilterbankRate = 44100.f / 512.f;
		float WindowSizeFreqHz = 200.f;
	} C;

	struct Parameters
	{
		float TonalTConstant = 0.1f;
		float TransientTConstant = 0.005f;
		float TonalThreshold = 2.f;
		float TransientThreshold = 1.25f;
		float PowerCompression = 0.3f;
	} P;

	struct Data
	{
		Eigen::ArrayXXf PowerMinSmooth, PowerFreqSmooth;
		float TonalLambda, TransientLambda;
		Eigen::ArrayXcf HFilter;
		void Reset()
		{
			PowerMinSmooth.setZero();
			PowerFreqSmooth.setZero();
			HFilter.setZero();
		}
		bool InitializeMemory(const Coefficients& c)
		{
			PowerMinSmooth.resize(c.NBands, c.NChannels);
			PowerFreqSmooth.resize(c.NBands, c.NChannels);
			HFilter.resize(c.NBands);
			return true;
		}
		size_t GetAllocatedMemorySize() const
		{
			return PowerMinSmooth.GetAllocatedMemorySize() + PowerFreqSmooth.GetAllocatedMemorySize() + HFilter.GetAllocatedMemorySize();
		}
		void OnParameterChange(const Parameters& p, const Coefficients& c)
		{
			TonalLambda = 1 - expf(-1.f / (c.FilterbankRate*p.TonalTConstant));
			TransientLambda = 1 - expf(-1.f / (c.FilterbankRate*p.TransientTConstant));
		}
	} D;

	DEFINEMEMBERALGORITHMS(2, FindMin, FFT);

	auto InitializeMembers()
	{
		auto cMin = FindMin.GetCoefficients();
		auto nFFT = (C.NBands - 1) * 2;
		auto length = C.WindowSizeFreqHz / C.SampleRate * nFFT;
		cMin.Length = static_cast<int>(std::round(length / 2.f)) * 2 + 1; // make odd
		cMin.NChannels = C.NChannels;
		bool flag = FindMin.Initialize(cMin);

		auto cFFT = FFT.GetCoefficients();
		cFFT.FFTSize = nFFT;
		flag &= FFT.Initialize(cFFT);
		
		// initialize D.HFilter using FFT
		Eigen::ArrayXf h = Eigen::ArrayXf::Zero(nFFT);
		h.head(length / 2) = 2.f / length; // averaging filter
		FFT.Process(h, D.HFilter);
		D.HFilter = D.HFilter.abs2(); // square averaging filter in freq domain

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
			// convolve with square averaging filter in freq domain
			Eigen::ArrayXf xPowerMirror(2 * (C.NBands - 1));
			xPowerMirror.head(C.NBands) = xPower.col(channel);
			xPowerMirror.tail(C.NBands - 2) = xPower.block(1, channel, C.NBands - 2, 1).colwise().reverse();
			Eigen::ArrayXcf convFilter(C.NBands);
			FFT.Process(xPowerMirror, convFilter);
			convFilter *= D.HFilter;
			FFT.Inverse(convFilter, xPowerMirror);
			Eigen::ArrayXf PowerFreqSmoothOld = D.PowerFreqSmooth.col(channel);
			D.PowerFreqSmooth.col(channel)  += D.TransientLambda * (xPowerMirror.head(C.NBands) - D.PowerFreqSmooth.col(channel));
			tonalDetection.col(channel) = D.PowerFreqSmooth.col(channel) < PowerFreqSmoothOld * P.TransientThreshold;
		}
		Eigen::ArrayXXf xMin(C.NBands, C.NChannels);
		FindMin.Process(xPower, xMin);
		xPower -= xMin + D.PowerMinSmooth;
		D.PowerMinSmooth += D.TonalLambda * xPower;
		tonalDetection = (tonalDetection && D.PowerMinSmooth > xMin * P.TonalThreshold);
	}

	void ProcessOff(Input xFreq, Output tonalDetection) { tonalDetection.setConstant(false); }
};

// #include "CriticalBands.h"
//
//// TODO: Support more than 1 channel
//class DetectTonalOld : public Base<DetectTonalOld, I::Real2D, O::Bool2D>
//{
//	friend Base<DetectTonalOld, I::Real2D, O::Bool2D>;
//
//public:
//	CriticalBands ConvertBands;
//	Eigen::ArrayXf GetTonalDifference() const { return D.Tonal; }
//
//private:
//	struct Coefficients 
//	{
//		int NBands = 513;
//		int NChannels = 2;
//		float SampleRate = 44100.f;
//		float FilterbankRate = 44100.f / 1024.f;
//		float FilterRelativeRange = .0227f; // Frequency range / Sample rate - Range to smooth transient difference over
//	} C;
//
//	struct Parameters 
//	{
//		float TonalTransientRatio = .26f; // Ratio of how much the tonal difference should be weighted against the transient difference. Can be calculated as:
//										  // MaxTimeDiff = sum(abs(win(1:NFFT-R)).^2) ./ sum(abs(win(1:NFFT)).^2);
//										  // EWIN = abs(fft(win, NFFT)). ^ 2;
//										  // MaxFreqDiff = EWIN(2) / EWIN(1);
//									      // TonalTransientRatio = MaxFreqDiff / MaxTimeDiff;
//		float TonalTC = 0.02f;
//		ConstrainedType<float> TonalThreshold = { 0.5f, 0.f, 1.f };
//	} P;
//
//	struct Data 
//	{
//		Eigen::ArrayXXf EOld;
//		Eigen::ArrayXXf Tonal;
//		float TonalLambda;
//		int TonalSize;
//		int NBandsCritical;
//		void Reset() 
//		{
//			EOld.setZero();
//			Tonal.setZero();
//		}
//		bool InitializeMemory(const Coefficients& c)
//		{ 
//			TonalSize = c.NBands - 1;
//			EOld.resize(c.NBands, c.NChannels);
//			Tonal.resize(TonalSize, c.NChannels);
//			return true;
//		}
//		size_t GetAllocatedMemorySize() const 
//		{ 
//			size_t size = Tonal.GetAllocatedMemorySize();
//			size += EOld.GetAllocatedMemorySize();
//			return  size;
//		}
//		void OnParameterChange(const Parameters& p, const Coefficients& c) 
//		{
//			TonalLambda = 1 - expf(-1.f / (c.FilterbankRate*p.TonalTC));
//		}
//	} D;
//
//	DEFINEMEMBERALGORITHMS(1, ConvertBands);
//
//	auto InitializeMembers()
//	{
//		auto c = ConvertBands.GetCoefficients();
//		c.NBands = C.NBands;
//		c.SampleRate = C.SampleRate;
//		bool flag = ConvertBands.Initialize(c);
//		D.NBandsCritical = ConvertBands.GetNBandsCritical();
//		return flag;
//	}
//
//	void ProcessOn(Input xEnergy, Output activityFlag) 
//	{
//		// freq difference
//		D.Tonal += D.TonalLambda * (xEnergy.bottomRows(D.TonalSize) - xEnergy.topRows(D.TonalSize) - D.Tonal);
//		Eigen::ArrayXXf tonalNum(C.NBands, C.NChannels);
//		tonalNum.topRows<1>() = D.Tonal.topRows<1>().abs();
//		tonalNum.bottomRows<1>() = D.Tonal.bottomRows<1>().abs();
//		tonalNum.middleRows(1, D.TonalSize - 1) = .5f*(D.Tonal.topRows(D.TonalSize - 1).abs() + D.Tonal.bottomRows(D.TonalSize - 1).abs());
//
//		// time difference
//		Eigen::ArrayXXf DiffCriticalBands(D.NBandsCritical, C.NChannels);
//		ConvertBands.Process((xEnergy - D.EOld).abs(), DiffCriticalBands);
//		Eigen::ArrayXXf tonalDen(C.NBands, C.NChannels);
//		ConvertBands.Inverse<float>(DiffCriticalBands, tonalDen);
//
//		// ratio
//		tonalNum *= P.TonalTransientRatio;
//		tonalNum /= (tonalDen + tonalNum); // Equals: y = tonalNum / tonalDen * P.TonalTransientRatio; yProbability = y/(y+1);
//
//		typedef Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> ArrayXXb;
//		activityFlag = (tonalNum > P.TonalThreshold).select(ArrayXXb::Constant(C.NBands, C.NChannels, true), ArrayXXb::Constant(C.NBands, C.NChannels, false));
//
//		D.EOld = xEnergy;
//	}
//
//	void ProcessOff(Input xFreq, Output activityFlag) { activityFlag.setConstant(false); }
//};