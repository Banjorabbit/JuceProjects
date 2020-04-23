#pragma once
#include "../BaseClasses/PureCRTP.h"
#include "CriticalBands.h"

// TODO: Support more than 1 channel
class DetectTonal : public Base<DetectTonal, I::Real2D, O::Bool2D>
{
	friend Base<DetectTonal, I::Real2D, O::Bool2D>;

public:
	CriticalBands ConvertBands;
	Eigen::ArrayXf GetTonalDifference() const { return D.Tonal; }

private:
	struct Coefficients 
	{
		int NBands = 513;
		int NChannels = 2;
		float SampleRate = 44100.f;
		float FilterbankRate = 44100.f / 1024.f;
		float FilterRelativeRange = .0227f; // Frequency range / Sample rate - Range to smooth transient difference over
	} C;

	struct Parameters 
	{
		float TonalTransientRatio = .26f; // Ratio of how much the tonal difference should be weighted against the transient difference. Can be calculated as:
										  // MaxTimeDiff = sum(abs(win(1:NFFT-R)).^2) ./ sum(abs(win(1:NFFT)).^2);
										  // EWIN = abs(fft(win, NFFT)). ^ 2;
										  // MaxFreqDiff = EWIN(2) / EWIN(1);
									      // TonalTransientRatio = MaxFreqDiff / MaxTimeDiff;
		float TonalTC = 0.02f;
		ConstrainedType<float> TonalThreshold = { 0.5f, 0.f, 1.f };
	} P;

	struct Data 
	{
		Eigen::ArrayXXf EOld;
		Eigen::ArrayXXf Tonal;
		float TonalLambda;
		int TonalSize;
		int NBandsCritical;
		void Reset() 
		{
			EOld.setZero();
			Tonal.setZero();
		}
		bool InitializeMemory(const Coefficients& c)
		{ 
			TonalSize = c.NBands - 1;
			EOld.resize(c.NBands, c.NChannels);
			Tonal.resize(TonalSize, c.NChannels);
			return true;
		}
		size_t GetAllocatedMemorySize() const 
		{ 
			size_t size = Tonal.GetAllocatedMemorySize();
			size += EOld.GetAllocatedMemorySize();
			return  size;
		}
		void OnParameterChange(const Parameters& p, const Coefficients& c) 
		{
			TonalLambda = 1 - expf(-1.f / (c.FilterbankRate*p.TonalTC));
		}
	} D;

	DEFINEMEMBERALGORITHMS(1, ConvertBands);

	auto InitializeMembers()
	{
		auto c = ConvertBands.GetCoefficients();
		c.NBands = C.NBands;
		c.SampleRate = C.SampleRate;
		bool flag = ConvertBands.Initialize(c);
		D.NBandsCritical = ConvertBands.GetNBandsCritical();
		return flag;
	}

	void ProcessOn(Input xEnergy, Output activityFlag) 
	{
		// freq difference
		D.Tonal += D.TonalLambda * (xEnergy.bottomRows(D.TonalSize) - xEnergy.topRows(D.TonalSize) - D.Tonal);
		Eigen::ArrayXXf tonalNum(C.NBands, C.NChannels);
		tonalNum.topRows<1>() = D.Tonal.topRows<1>().abs();
		tonalNum.bottomRows<1>() = D.Tonal.bottomRows<1>().abs();
		tonalNum.middleRows(1, D.TonalSize - 1) = .5f*(D.Tonal.topRows(D.TonalSize - 1).abs() + D.Tonal.bottomRows(D.TonalSize - 1).abs());

		// time difference
		Eigen::ArrayXXf DiffCriticalBands(D.NBandsCritical, C.NChannels);
		ConvertBands.Process((xEnergy - D.EOld).abs(), DiffCriticalBands);
		Eigen::ArrayXXf tonalDen(C.NBands, C.NChannels);
		ConvertBands.Inverse<float>(DiffCriticalBands, tonalDen);

		// ratio
		tonalNum *= P.TonalTransientRatio;
		tonalNum /= (tonalDen + tonalNum); // Equals: y = tonalNum / tonalDen * P.TonalTransientRatio; yProbability = y/(y+1);

		typedef Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> ArrayXXb;
		activityFlag = (tonalNum > P.TonalThreshold).select(ArrayXXb::Constant(C.NBands, C.NChannels, true), ArrayXXb::Constant(C.NBands, C.NChannels, false));

		D.EOld = xEnergy;
	}

	void ProcessOff(Input xFreq, Output activityFlag) { activityFlag.setConstant(false); }
};