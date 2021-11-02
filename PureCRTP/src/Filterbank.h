#pragma once
#include "BaseClasses/PureCRTP.h"
#include "FFT.h"

class FilterbankAnalysis : public Base<FilterbankAnalysis, I::Real2D, O::Complex2D>
{
	friend Base<FilterbankAnalysis, I::Real2D, O::Complex2D>;

public:
	FFTReal FFT;
	void SetWindow(const Eigen::Ref<const Eigen::ArrayXf>& window) 
	{
		if ((window.size() == C.FrameSize) & (P.WindowType == P.UserDefined)) { D.Window = window; }
	}

	Eigen::ArrayXf GetWindow() const { return D.Window; }

	static Eigen::ArrayXf GetHannWindow(const int size) { return .5f * (1.f - Eigen::ArrayXf::LinSpaced(size, 0, 2.f * static_cast<float>(PI) * (size - 1) / size).cos()); };

private:
	struct Coefficients
	{
		int NChannels = 2;
		int FrameSize = 512;
		int FFTSize = 512;
		int BufferSize = 128;
	} C;

	struct Parameters
	{
		enum WindowTypes { HannWindow , SqrtHannWindow, Rectangular, UserDefined};
		WindowTypes WindowType = HannWindow;
		float Gain = 1;
	} P;

	struct Data
	{
		Eigen::ArrayXf Window, FFTBuffer; // advantage of FFTBuffer being allocated on heap is that zeropadding is kept between calls
		Eigen::ArrayXXf TimeInputBuffer;
		int Overlap, NFolds, MaxSize;
		void Reset()
		{
			TimeInputBuffer.setZero();
			FFTBuffer.setZero();
		}
		bool InitializeMemory(const Coefficients& c)
		{
			Overlap = c.FrameSize - c.BufferSize;
			NFolds = static_cast<int>(std::ceil(static_cast<float>(c.FrameSize) / c.FFTSize));
			MaxSize = c.FFTSize * NFolds;
			Window.resize(c.FrameSize);
			FFTBuffer.resize(MaxSize);
			TimeInputBuffer.resize(c.FrameSize, c.NChannels);
			return true;
		}
		size_t GetAllocatedMemorySize() const
		{
			size_t size = Window.GetAllocatedMemorySize();
			size += TimeInputBuffer.GetAllocatedMemorySize();
			size += FFTBuffer.GetAllocatedMemorySize();
			return size;
		}
		void OnParameterChange(const Parameters& p, const Coefficients& c)
		{
			switch (p.WindowType)
			{
			case Parameters::HannWindow:
				Window = GetHannWindow(c.FrameSize);
				break;
			case Parameters::SqrtHannWindow:
				Window = GetHannWindow(c.FrameSize).sqrt();
				break;
			case Parameters::Rectangular:
				Window.setOnes();
				break;
			case Parameters::UserDefined:
				break;
			default:
				Window.setZero();
			}
			Window *= p.Gain;
		}
	} D;

	DEFINEMEMBERALGORITHMS(1, FFT);

	auto InitializeMembers()
	{
		auto sFFT = FFT.GetSetup();
		sFFT.Coefficients.FFTSize = C.FFTSize;
		return FFT.Initialize(sFFT);
	}

	void ProcessOn(Input xTime, Output yFreq)
	{
		D.TimeInputBuffer.topRows(D.Overlap) = D.TimeInputBuffer.bottomRows(D.Overlap);
		D.TimeInputBuffer.bottomRows(C.BufferSize) = xTime;
		for (auto channel = 0; channel < xTime.cols();channel++)
		{
			D.FFTBuffer.head(C.FrameSize) = D.TimeInputBuffer.col(channel) * D.Window;
			for (auto j = 1; j < D.NFolds; j++)
			{
				D.FFTBuffer.head(C.FFTSize) += D.FFTBuffer.segment(j*C.FFTSize, C.FFTSize);
			}
			FFT.Process(D.FFTBuffer.head(C.FFTSize), yFreq.col(channel));
		}
	}

	void ProcessOff(Input xTime, Output yFreq) { yFreq.setZero(); }
};

class FilterbankSynthesis : public Base<FilterbankSynthesis, I::Complex2D, O::Real2D>
{
	friend Base<FilterbankSynthesis, I::Complex2D, O::Real2D>;

public:
	FFTRealInverse FFTInv;
	Eigen::ArrayXf GetWindow() const { return D.Window; }

	static Eigen::ArrayXf GetHannWindow(const int size) { return .5f * (1.f - Eigen::ArrayXf::LinSpaced(size, 0, 2.f * static_cast<float>(PI) * (size - 1) / size).cos()); };

private:
	struct Coefficients
	{
		int NChannels = 2;
		int FrameSize = 512;
		int FFTSize = 512;
		int BufferSize = 128;
	} C;

	struct Parameters
	{
		enum WindowTypes { HannWindow , SqrtHannWindow, Rectangular};
		WindowTypes WindowType = HannWindow;
		float Gain = 1;
	} P;

	struct Data
	{
		Eigen::ArrayXf Window, FftBuffer;
		Eigen::ArrayXXf TimeOutputBuffer;
		int Overlap, NFolds, MaxSize;
		void Reset()
		{
			FftBuffer.setZero();
			TimeOutputBuffer.setZero();
		}
		bool InitializeMemory(const Coefficients& c)
		{
			Overlap = c.FrameSize - c.BufferSize;
			NFolds = static_cast<int>(std::ceil(static_cast<float>(c.FrameSize) / c.FFTSize));
			MaxSize = c.FFTSize * NFolds;
			FftBuffer.resize(MaxSize);
			Window.resize(c.FrameSize);
			TimeOutputBuffer.resize(c.FrameSize, c.NChannels);
			return true;
		}
		size_t GetAllocatedMemorySize() const
		{
			size_t size = Window.GetAllocatedMemorySize();
			size += FftBuffer.GetAllocatedMemorySize();
			size += TimeOutputBuffer.GetAllocatedMemorySize();
			return size;
		}
		void OnParameterChange(const Parameters& p, const Coefficients& c)
		{
			switch (p.WindowType)
			{
			case Parameters::HannWindow:
				Window = GetHannWindow(c.FrameSize);
				break;
			case Parameters::SqrtHannWindow:
				Window = GetHannWindow(c.FrameSize).sqrt();
				break;
			case Parameters::Rectangular:
				Window.setOnes();
				break;
			default:
				Window.setZero();
			}
			Window *= p.Gain;
		}
	} D;

	DEFINEMEMBERALGORITHMS(1, FFTInv);

	auto InitializeMembers()
	{
		auto sFFT = FFTInv.GetSetup();
		sFFT.Coefficients.FFTSize = C.FFTSize;
		return FFTInv.Initialize(sFFT);
	}

	void ProcessOn(Input xFreq, Output yTime)
	{
		for (auto channel = 0; channel < xFreq.cols();channel++)
		{
			FFTInv.Process(xFreq.col(channel), D.FftBuffer.head(C.FFTSize));
			for (auto j = 1; j < D.NFolds;j++)
			{
				D.FftBuffer.segment(j*C.FFTSize, C.FFTSize) = D.FftBuffer.head(C.FFTSize);
			}
			D.TimeOutputBuffer.col(channel) += (D.FftBuffer.head(C.FrameSize) * D.Window);
		}
		yTime = D.TimeOutputBuffer.topRows(C.BufferSize);
		D.TimeOutputBuffer.topRows(D.Overlap) = D.TimeOutputBuffer.bottomRows(D.Overlap);
		D.TimeOutputBuffer.bottomRows(C.BufferSize) = 0.f;
	}

	void ProcessOff(Input xFreq, Output yTime) { yTime.setZero(); }
};