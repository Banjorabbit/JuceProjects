#pragma once
#include "../BaseClasses/PureCRTP.h"

class TonalDetection : public Base<TonalDetection, I::Real, O::Real>
{
	friend Base<TonalDetection, I::Real, O::Real>;

public:
	Eigen::ArrayXf GetTonalDifference() const { return D.Tonal; }

private:
	struct Coefficients 
	{
		int NBands = 513;
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
	} P;

	struct Data 
	{
		Eigen::ArrayXf EOld;
		Eigen::ArrayXf Tonal;
		Eigen::ArrayXf FilterTransient;
		float TonalLambda;
		int TonalSize, TransientSize, FilterOrder, FilterHalfOrder;
		void Reset() 
		{
			EOld.setZero();
			Tonal.setZero();
		}
		bool InitializeMemory(const Coefficients& c)
		{
			FilterHalfOrder = static_cast<int>(c.FilterRelativeRange*(c.NBands - 1.f));
			FilterOrder = FilterHalfOrder*2 + 1; 
			TonalSize = c.NBands - 1;
			Eigen::ArrayXf filterHalf = Eigen::ArrayXf::LinSpaced(FilterHalfOrder+1, 1.f, FilterHalfOrder + 1.f);
			FilterTransient.resize(FilterOrder);
			FilterTransient.head(FilterHalfOrder+1) = filterHalf;
			FilterTransient.tail(FilterHalfOrder) = filterHalf.head(FilterHalfOrder).reverse();
			FilterTransient /= FilterTransient.sum();
			EOld.resize(c.NBands);
			Tonal.resize(TonalSize);
			return true;
		}
		size_t GetAllocatedMemorySize() const 
		{ 
			size_t size = Tonal.GetAllocatedMemorySize();
			size += EOld.GetAllocatedMemorySize();
			size += FilterTransient.GetAllocatedMemorySize();
			return  size;
		}
		void OnParameterChange(const Parameters& p, const Coefficients& c) 
		{
			TonalLambda = 1 - expf(-1.f / (c.FilterbankRate*p.TonalTC));
		}
	} D;

	void ProcessOn(Input xEnergy, Output yMag) 
	{
		const Eigen::ArrayXf tonalDifference = (xEnergy.tail(D.TonalSize) - xEnergy.head(D.TonalSize));
		D.Tonal += D.TonalLambda * (tonalDifference - D.Tonal);
		Eigen::ArrayXf transientDifference(C.NBands + 2 * D.FilterHalfOrder);
		transientDifference.segment(D.FilterHalfOrder, C.NBands) = (xEnergy - D.EOld).abs();
		transientDifference.head(D.FilterHalfOrder) = transientDifference.segment(D.FilterHalfOrder + 1, D.FilterHalfOrder).reverse();
		transientDifference.tail(D.FilterHalfOrder) = transientDifference.segment(D.TonalSize, D.FilterHalfOrder).reverse();	
		Eigen::ArrayXf tonalDen = Eigen::ArrayXf::Constant(C.NBands, 1e-20f);
		for (auto i = 0; i < D.FilterOrder; i++)
		{
			tonalDen += transientDifference.segment(i, C.NBands) * D.FilterTransient(i);
		}
		Eigen::ArrayXf tonalNum(C.NBands);
		tonalNum(0) = std::abs(D.Tonal(0));
		tonalNum(D.TonalSize) = std::abs(D.Tonal(D.TonalSize - 1));
		tonalNum.segment(1, D.TonalSize - 1) = .5f*(D.Tonal.head(D.TonalSize - 1).abs() + D.Tonal.tail(D.TonalSize - 1).abs());
		yMag = tonalNum / tonalDen * P.TonalTransientRatio;
		D.EOld = xEnergy;
	}

	void ProcessOff(Input xFreq, Output yMag) { yMag.setZero(); }
};