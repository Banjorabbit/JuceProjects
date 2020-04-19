#pragma once
#include "../BaseClasses/PureCRTP.h"

struct O::IIR2ndTimeVaryingFilter
{
	Real2D HighPass;
	Real2D LowPass;
	Real2D BandPass;
};

class IIR2ndTimeVaryingFilter : public Base<IIR2ndTimeVaryingFilter, I::Real2D, O::IIR2ndTimeVaryingFilter>
{
	friend Base<IIR2ndTimeVaryingFilter, I::Real2D, O::IIR2ndTimeVaryingFilter>;

	struct Coefficients
	{
		float SampleRate = 16e3f;
		int NChannels = 2;
	} C;

	struct Parameters
	{
		ConstrainedType<float> Cutoff = { 1000.f, 1.f, 20e3f };
		ConstrainedType<float> Q = { static_cast<float>(BUTTERWORTH_Q), 1e-3f, 100.f };
	} P;

	struct Data
	{
		Eigen::ArrayXf State1, State2;
		float c0,c1,c2,c3;
		void Reset()
		{
			State1.setZero();
			State2.setZero();
		}
		bool InitializeMemory(const Coefficients& c)
		{
			State1.resize(c.NChannels);
			State2.resize(c.NChannels);
			return true;
		}
		size_t GetAllocatedMemorySize() const
		{
			auto size = State1.GetAllocatedMemorySize();
			size += State2.GetAllocatedMemorySize();
			return size;
		}
		void OnParameterChange(const Parameters& p, const Coefficients& c)
		{
			const auto g = static_cast<float>(std::tan(PI*p.Cutoff / c.SampleRate)); 
			c1 = g / p.Q;
			c2 = g * p.Q;
			c3 = c2 + 1.f;
			c0 = 1.f / (1.f + c1 * (c2 + 1.f));
		}
	} D;

	void ProcessOn(Input xTime, Output yTime) 
	{
		for (auto sample = 0; sample < xTime.rows(); sample++)
		{
			for (auto channel = 0; channel < C.NChannels; channel++)
			{
				yTime.HighPass(sample, channel) = D.c0 * (xTime(sample, channel) - D.State2(channel) - D.c3 * D.State1(channel));
				const auto x1 = D.c1 * yTime.HighPass(sample, channel);
				yTime.BandPass(sample, channel) = x1 + D.State1(channel);
				const auto x2 = D.c2 * yTime.BandPass(sample, channel);
				yTime.LowPass(sample, channel) = x2 + D.State2(channel);

				D.State1(channel) = x1 + yTime.BandPass(sample, channel);
				D.State2(channel) = x2 + yTime.LowPass(sample, channel);
			}
			
		}
		// combine to make Notch, Bell, HighShelf, LowShelf:
		// D.g0 = 10.^(gain./20);
		// D.g1 = D.g0 - 1;
		//yTime.Notch = xTime - yTime.BandPass;
		//yTime.Bell = xTime + D.g1 * yTime.BandPass;
		//yTime.HighShelf = D.g0 * (yTime.BandPass + yTime.HighPass) + yTime.LowPass;
		//yTime.LowShelf = D.g0 * (yTime.BandPass + yTime.LowPass) + yTime.HighPass;
	}

	void ProcessOff(Input xTime, Output yTime) 
	{ 
		yTime.HighPass = xTime;
		yTime.LowPass = xTime;
		yTime.BandPass = xTime;
	}
};

class IIR2ndTimeVaryingLowPassFilter : public Base<IIR2ndTimeVaryingLowPassFilter>
{
	friend Base<IIR2ndTimeVaryingLowPassFilter>;
public:
	IIR2ndTimeVaryingFilter Filter;
private:
	struct Coefficients	{ float SampleRate = 16e3f;	int NChannels = 2; } C;
	struct Parameters {} P;
	struct Data 
	{
		void Reset() {}
		bool InitializeMemory(const Coefficients& c) { return true; }
		size_t GetAllocatedMemorySize() const { return 0; }
		void OnParameterChange(const Parameters& p, const Coefficients& c) {}
	} D;
	DEFINEMEMBERALGORITHMS(1, Filter);
	auto InitializeMembers() { auto c = Filter.GetCoefficients(); c.NChannels = C.NChannels; c.SampleRate = C.SampleRate; return Filter.Initialize(c); }
	void ProcessOn(Input xTime, Output yTime) { Eigen::ArrayXXf bandPass(xTime.rows(), xTime.cols()), highPass(xTime.rows(), xTime.cols()); Filter.Process(xTime, { highPass, yTime, bandPass }); }
	void ProcessOff(Input xTime, Output yTime) { yTime = xTime; }
};

class IIR2ndTimeVaryingHighPassFilter : public Base<IIR2ndTimeVaryingHighPassFilter>
{
	friend Base<IIR2ndTimeVaryingHighPassFilter>;
public:
	IIR2ndTimeVaryingFilter Filter;
private:
	struct Coefficients { float SampleRate = 16e3f;	int NChannels = 2; } C;
	struct Parameters {} P;
	struct Data
	{
		void Reset() {}
		bool InitializeMemory(const Coefficients& c) { return true; }
		size_t GetAllocatedMemorySize() const { return 0; }
		void OnParameterChange(const Parameters& p, const Coefficients& c) {}
	} D;
	DEFINEMEMBERALGORITHMS(1, Filter);
	auto InitializeMembers() { auto c = Filter.GetCoefficients(); c.NChannels = C.NChannels; c.SampleRate = C.SampleRate; return Filter.Initialize(c); }
	void ProcessOn(Input xTime, Output yTime) { Eigen::ArrayXXf bandPass(xTime.rows(), xTime.cols()), lowPass(xTime.rows(), xTime.cols()); Filter.Process(xTime, { yTime, lowPass, bandPass }); }
	void ProcessOff(Input xTime, Output yTime) { yTime = xTime; }
};

class IIR2ndTimeVaryingBandPassFilter : public Base<IIR2ndTimeVaryingBandPassFilter>
{
	friend Base<IIR2ndTimeVaryingBandPassFilter>;
public:
	IIR2ndTimeVaryingFilter Filter;
private:
	struct Coefficients { float SampleRate = 16e3f;	int NChannels = 2; } C;
	struct Parameters {} P;
	struct Data
	{
		void Reset() {}
		bool InitializeMemory(const Coefficients& c) { return true; }
		size_t GetAllocatedMemorySize() const { return 0; }
		void OnParameterChange(const Parameters& p, const Coefficients& c) {}
	} D;
	DEFINEMEMBERALGORITHMS(1, Filter);
	auto InitializeMembers() { auto c = Filter.GetCoefficients(); c.NChannels = C.NChannels; c.SampleRate = C.SampleRate; return Filter.Initialize(c); }
	void ProcessOn(Input xTime, Output yTime) { Eigen::ArrayXXf highPass(xTime.rows(), xTime.cols()), lowPass(xTime.rows(), xTime.cols()); Filter.Process(xTime, { highPass, lowPass, yTime }); }
	void ProcessOff(Input xTime, Output yTime) { yTime = xTime; }
};

class IIR2ndTimeVaryingNotchFilter : public Base<IIR2ndTimeVaryingNotchFilter>
{
	friend Base<IIR2ndTimeVaryingNotchFilter>;
public:
	IIR2ndTimeVaryingFilter Filter;
private:
	struct Coefficients { float SampleRate = 16e3f;	int NChannels = 2; } C;
	struct Parameters {} P;
	struct Data
	{
		void Reset() {}
		bool InitializeMemory(const Coefficients& c) { return true; }
		size_t GetAllocatedMemorySize() const { return 0; }
		void OnParameterChange(const Parameters& p, const Coefficients& c) {}
	} D;
	DEFINEMEMBERALGORITHMS(1, Filter);
	auto InitializeMembers() { auto c = Filter.GetCoefficients(); c.NChannels = C.NChannels; c.SampleRate = C.SampleRate; return Filter.Initialize(c); }
	void ProcessOn(Input xTime, Output yTime) 
	{ 
		Eigen::ArrayXXf highPass(xTime.rows(), xTime.cols()), lowPass(xTime.rows(), xTime.cols()), bandPass(xTime.rows(), xTime.cols()); 
		Filter.Process(xTime, { highPass, lowPass, bandPass }); 
		yTime = xTime - bandPass;
	}
	void ProcessOff(Input xTime, Output yTime) { yTime = xTime; }
};

class IIR2ndTimeVaryingBellFilter : public Base<IIR2ndTimeVaryingBellFilter>
{
	friend Base<IIR2ndTimeVaryingBellFilter>;
public:
	IIR2ndTimeVaryingFilter Filter;
private:
	struct Coefficients { float SampleRate = 16e3f;	int NChannels = 2; } C;
	struct Parameters { float GaindB = 0.f; } P;
	struct Data
	{
		float Gain;
		void Reset() {}
		bool InitializeMemory(const Coefficients& c) { return true; }
		size_t GetAllocatedMemorySize() const { return 0; }
		void OnParameterChange(const Parameters& p, const Coefficients& c) { Gain = powf(10.f, p.GaindB * .05f) - 1.f; }
	} D;
	DEFINEMEMBERALGORITHMS(1, Filter);
	auto InitializeMembers() { auto c = Filter.GetCoefficients(); c.NChannels = C.NChannels; c.SampleRate = C.SampleRate; return Filter.Initialize(c); }
	void ProcessOn(Input xTime, Output yTime)
	{
		Eigen::ArrayXXf highPass(xTime.rows(), xTime.cols()), lowPass(xTime.rows(), xTime.cols()), bandPass(xTime.rows(), xTime.cols());
		Filter.Process(xTime, { highPass, lowPass, bandPass });
		yTime = xTime + D.Gain * bandPass;
	}
	void ProcessOff(Input xTime, Output yTime) { yTime = xTime; }
};

class IIR2ndTimeVaryingHighShelfFilter : public Base<IIR2ndTimeVaryingHighShelfFilter>
{
	friend Base<IIR2ndTimeVaryingHighShelfFilter>;
public:
	IIR2ndTimeVaryingFilter Filter;
private:
	struct Coefficients { float SampleRate = 16e3f;	int NChannels = 2; } C;
	struct Parameters { float GaindB = 0.f; } P;
	struct Data
	{
		float Gain;
		void Reset() {}
		bool InitializeMemory(const Coefficients& c) { return true; }
		size_t GetAllocatedMemorySize() const { return 0; }
		void OnParameterChange(const Parameters& p, const Coefficients& c) { Gain = powf(10.f, p.GaindB * .05f); }
	} D;
	DEFINEMEMBERALGORITHMS(1, Filter);
	auto InitializeMembers() { auto c = Filter.GetCoefficients(); c.NChannels = C.NChannels; c.SampleRate = C.SampleRate; return Filter.Initialize(c); }
	void ProcessOn(Input xTime, Output yTime)
	{
		Eigen::ArrayXXf highPass(xTime.rows(), xTime.cols()), lowPass(xTime.rows(), xTime.cols()), bandPass(xTime.rows(), xTime.cols());
		Filter.Process(xTime, { highPass, lowPass, bandPass });
		yTime = D.Gain * (bandPass + highPass) + lowPass;
	}
	void ProcessOff(Input xTime, Output yTime) { yTime = xTime; }
};

class IIR2ndTimeVaryingLowShelfFilter : public Base<IIR2ndTimeVaryingLowShelfFilter>
{
	friend Base<IIR2ndTimeVaryingLowShelfFilter>;
public:
	IIR2ndTimeVaryingFilter Filter;
private:
	struct Coefficients { float SampleRate = 16e3f;	int NChannels = 2; } C;
	struct Parameters { float GaindB = 0.f; } P;
	struct Data
	{
		float Gain;
		void Reset() {}
		bool InitializeMemory(const Coefficients& c) { return true; }
		size_t GetAllocatedMemorySize() const { return 0; }
		void OnParameterChange(const Parameters& p, const Coefficients& c) { Gain = powf(10.f, p.GaindB * .05f); }
	} D;
	DEFINEMEMBERALGORITHMS(1, Filter);
	auto InitializeMembers() { auto c = Filter.GetCoefficients(); c.NChannels = C.NChannels; c.SampleRate = C.SampleRate; return Filter.Initialize(c); }
	void ProcessOn(Input xTime, Output yTime)
	{
		Eigen::ArrayXXf highPass(xTime.rows(), xTime.cols()), lowPass(xTime.rows(), xTime.cols()), bandPass(xTime.rows(), xTime.cols());
		Filter.Process(xTime, { highPass, lowPass, bandPass });
		yTime = D.Gain * (bandPass + lowPass) + highPass;
	}
	void ProcessOff(Input xTime, Output yTime) { yTime = xTime; }
};

// ------------------ Streaming classes ------------------------------------------

class IIR2ndTimeVaryingLowPassFilterStreaming : public AsynchronousStreaming<IIR2ndTimeVaryingLowPassFilterStreaming, IIR2ndTimeVaryingLowPassFilter>
{
	friend AsynchronousStreaming<IIR2ndTimeVaryingLowPassFilterStreaming, IIR2ndTimeVaryingLowPassFilter>;

	int GetLatencySamples(const decltype(Algo.GetCoefficients())& c) const { return 0; }
	int GetNChannelsOut(const int nChannels) const { return nChannels; }

	int UpdateCoefficients(decltype(Algo.GetCoefficients())& c, const float sampleRate, const int nChannels, const int bufferSizeExpected, std::vector<int> bufferSizesSuggested) const
	{
		c.NChannels = nChannels;
		c.SampleRate = sampleRate;
		return bufferSizeExpected;
	}
};

class IIR2ndTimeVaryingHighPassFilterStreaming : public AsynchronousStreaming<IIR2ndTimeVaryingHighPassFilterStreaming, IIR2ndTimeVaryingHighPassFilter>
{
	friend AsynchronousStreaming<IIR2ndTimeVaryingHighPassFilterStreaming, IIR2ndTimeVaryingHighPassFilter>;

	int GetLatencySamples(const decltype(Algo.GetCoefficients())& c) const { return 0; }
	int GetNChannelsOut(const int nChannels) const { return nChannels; }

	int UpdateCoefficients(decltype(Algo.GetCoefficients())& c, const float sampleRate, const int nChannels, const int bufferSizeExpected, std::vector<int> bufferSizesSuggested) const
	{
		c.NChannels = nChannels;
		c.SampleRate = sampleRate;
		return bufferSizeExpected;
	}
};

class IIR2ndTimeVaryingBandPassFilterStreaming : public AsynchronousStreaming<IIR2ndTimeVaryingBandPassFilterStreaming, IIR2ndTimeVaryingBandPassFilter>
{
	friend AsynchronousStreaming<IIR2ndTimeVaryingBandPassFilterStreaming, IIR2ndTimeVaryingBandPassFilter>;

	int GetLatencySamples(const decltype(Algo.GetCoefficients())& c) const { return 0; }
	int GetNChannelsOut(const int nChannels) const { return nChannels; }

	int UpdateCoefficients(decltype(Algo.GetCoefficients())& c, const float sampleRate, const int nChannels, const int bufferSizeExpected, std::vector<int> bufferSizesSuggested) const
	{
		c.NChannels = nChannels;
		c.SampleRate = sampleRate;
		return bufferSizeExpected;
	}
};

class IIR2ndTimeVaryingNotchFilterStreaming : public AsynchronousStreaming<IIR2ndTimeVaryingNotchFilterStreaming, IIR2ndTimeVaryingNotchFilter>
{
	friend AsynchronousStreaming<IIR2ndTimeVaryingNotchFilterStreaming, IIR2ndTimeVaryingNotchFilter>;

	int GetLatencySamples(const decltype(Algo.GetCoefficients())& c) const { return 0; }
	int GetNChannelsOut(const int nChannels) const { return nChannels; }

	int UpdateCoefficients(decltype(Algo.GetCoefficients())& c, const float sampleRate, const int nChannels, const int bufferSizeExpected, std::vector<int> bufferSizesSuggested) const
	{
		c.NChannels = nChannels;
		c.SampleRate = sampleRate;
		return bufferSizeExpected;
	}
};

class IIR2ndTimeVaryingBellFilterStreaming : public AsynchronousStreaming<IIR2ndTimeVaryingBellFilterStreaming, IIR2ndTimeVaryingBellFilter>
{
	friend AsynchronousStreaming<IIR2ndTimeVaryingBellFilterStreaming, IIR2ndTimeVaryingBellFilter>;

	int GetLatencySamples(const decltype(Algo.GetCoefficients())& c) const { return 0; }
	int GetNChannelsOut(const int nChannels) const { return nChannels; }

	int UpdateCoefficients(decltype(Algo.GetCoefficients())& c, const float sampleRate, const int nChannels, const int bufferSizeExpected, std::vector<int> bufferSizesSuggested) const
	{
		c.NChannels = nChannels;
		c.SampleRate = sampleRate;
		return bufferSizeExpected;
	}
};

class IIR2ndTimeVaryingHighShelfFilterStreaming : public AsynchronousStreaming<IIR2ndTimeVaryingHighShelfFilterStreaming, IIR2ndTimeVaryingHighShelfFilter>
{
	friend AsynchronousStreaming<IIR2ndTimeVaryingHighShelfFilterStreaming, IIR2ndTimeVaryingHighShelfFilter>;

	int GetLatencySamples(const decltype(Algo.GetCoefficients())& c) const { return 0; }
	int GetNChannelsOut(const int nChannels) const { return nChannels; }

	int UpdateCoefficients(decltype(Algo.GetCoefficients())& c, const float sampleRate, const int nChannels, const int bufferSizeExpected, std::vector<int> bufferSizesSuggested) const
	{
		c.NChannels = nChannels;
		c.SampleRate = sampleRate;
		return bufferSizeExpected;
	}
};

class IIR2ndTimeVaryingLowShelfFilterStreaming : public AsynchronousStreaming<IIR2ndTimeVaryingLowShelfFilterStreaming, IIR2ndTimeVaryingLowShelfFilter>
{
	friend AsynchronousStreaming<IIR2ndTimeVaryingLowShelfFilterStreaming, IIR2ndTimeVaryingLowShelfFilter>;

	int GetLatencySamples(const decltype(Algo.GetCoefficients())& c) const { return 0; }
	int GetNChannelsOut(const int nChannels) const { return nChannels; }

	int UpdateCoefficients(decltype(Algo.GetCoefficients())& c, const float sampleRate, const int nChannels, const int bufferSizeExpected, std::vector<int> bufferSizesSuggested) const
	{
		c.NChannels = nChannels;
		c.SampleRate = sampleRate;
		return bufferSizeExpected;
	}
};