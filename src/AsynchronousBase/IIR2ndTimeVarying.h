#pragma once
#include "../BaseClasses/PureCRTP.h"
#include "../StateVariableFilter.h"

class IIR2ndTimeVaryingLowPassFilter : public AsynchronousBase<IIR2ndTimeVaryingLowPassFilter>
{
	friend Base<IIR2ndTimeVaryingLowPassFilter>;
public:
	StateVariableFilter Filter;
private:
	struct Coefficients 
	{ 
		float SampleRate = 16e3f;
		int NChannelsIn = 2; 
		int BufferSize = 128;
		AsynchronousBufferType AsynchronousBuffer = VARIABLE_SIZE;
	} C;

	struct Parameters {} P;

	struct Data 
	{
		void Reset() {}
		bool InitializeMemory(const Coefficients& c) { return true; }
		size_t GetAllocatedMemorySize() const { return 0; }
		void OnParameterChange(const Parameters& p, const Coefficients& c) {}
	} D;

	DEFINEMEMBERALGORITHMS(1, Filter);

	auto InitializeMembers() 
	{ 
		auto c = Filter.GetCoefficients(); 
		c.NChannels = C.NChannelsIn; 
		c.BufferSize = C.BufferSize; 
		c.SampleRate = C.SampleRate; 
		return Filter.Initialize(c); 
	}

	void ProcessOn(Input xTime, Output yTime) 
	{ 
		Eigen::ArrayXXf bandPass(xTime.rows(), xTime.cols()), highPass(xTime.rows(), xTime.cols()); 
		Filter.Process(xTime, { highPass, yTime, bandPass }); 
	}

	void ProcessOff(Input xTime, Output yTime) { yTime = xTime; }
};

class IIR2ndTimeVaryingHighPassFilter : public AsynchronousBase<IIR2ndTimeVaryingHighPassFilter>
{
	friend Base<IIR2ndTimeVaryingHighPassFilter>;
public:
	StateVariableFilter Filter;
private:
	struct Coefficients 
	{ 
		float SampleRate = 16e3f;
		int BufferSize = 128; 
		int NChannelsIn = 2; 
		AsynchronousBufferType AsynchronousBuffer = VARIABLE_SIZE;
	} C;

	struct Parameters {} P;

	struct Data
	{
		void Reset() {}
		bool InitializeMemory(const Coefficients& c) { return true; }
		size_t GetAllocatedMemorySize() const { return 0; }
		void OnParameterChange(const Parameters& p, const Coefficients& c) {}
	} D;

	DEFINEMEMBERALGORITHMS(1, Filter);

	auto InitializeMembers() 
	{ 
		auto c = Filter.GetCoefficients(); 
		c.NChannels = C.NChannelsIn; 
		c.BufferSize = C.BufferSize; 
		c.SampleRate = C.SampleRate; 
		return Filter.Initialize(c); 
	}

	void ProcessOn(Input xTime, Output yTime) 
	{ 
		Eigen::ArrayXXf bandPass(xTime.rows(), xTime.cols()), lowPass(xTime.rows(), xTime.cols());
		Filter.Process(xTime, { yTime, lowPass, bandPass }); 
	}

	void ProcessOff(Input xTime, Output yTime) { yTime = xTime; }
};

class IIR2ndTimeVaryingBandPassFilter : public AsynchronousBase<IIR2ndTimeVaryingBandPassFilter>
{
	friend Base<IIR2ndTimeVaryingBandPassFilter>;
public:
	StateVariableFilter Filter;
private:
	struct Coefficients 
	{ 
		float SampleRate = 16e3f;
		int NChannelsIn = 2;
		int BufferSize = 128;
		AsynchronousBufferType AsynchronousBuffer = VARIABLE_SIZE;
	} C;

	struct Parameters {} P;

	struct Data
	{
		void Reset() {}
		bool InitializeMemory(const Coefficients& c) { return true; }
		size_t GetAllocatedMemorySize() const { return 0; }
		void OnParameterChange(const Parameters& p, const Coefficients& c) {}
	} D;

	DEFINEMEMBERALGORITHMS(1, Filter);

	auto InitializeMembers() 
	{
		auto c = Filter.GetCoefficients(); 
		c.NChannels = C.NChannelsIn; 
		c.BufferSize = C.BufferSize; 
		c.SampleRate = C.SampleRate; 
		return Filter.Initialize(c);
	}

	void ProcessOn(Input xTime, Output yTime) 
	{ 
		Eigen::ArrayXXf highPass(xTime.rows(), xTime.cols()), lowPass(xTime.rows(), xTime.cols());
		Filter.Process(xTime, { highPass, lowPass, yTime }); 
	}

	void ProcessOff(Input xTime, Output yTime) { yTime = xTime; }
};

class IIR2ndTimeVaryingNotchFilter : public AsynchronousBase<IIR2ndTimeVaryingNotchFilter>
{
	friend Base<IIR2ndTimeVaryingNotchFilter>;
public:
	StateVariableFilter Filter;
private:
	struct Coefficients
	{
		float SampleRate = 16e3f;
		int NChannelsIn = 2;
		int BufferSize = 128;
		AsynchronousBufferType AsynchronousBuffer = VARIABLE_SIZE;
	} C;

	struct Parameters {} P;

	struct Data
	{
		void Reset() {}
		bool InitializeMemory(const Coefficients& c) { return true; }
		size_t GetAllocatedMemorySize() const { return 0; }
		void OnParameterChange(const Parameters& p, const Coefficients& c) {}
	} D;

	DEFINEMEMBERALGORITHMS(1, Filter);

	auto InitializeMembers() 
	{
		auto c = Filter.GetCoefficients();
		c.NChannels = C.NChannelsIn;
		c.BufferSize = C.BufferSize; 
		c.SampleRate = C.SampleRate; 
		return Filter.Initialize(c);
	}

	void ProcessOn(Input xTime, Output yTime) 
	{ 
		Eigen::ArrayXXf highPass(xTime.rows(), xTime.cols()), lowPass(xTime.rows(), xTime.cols()), bandPass(xTime.rows(), xTime.cols()); 
		Filter.Process(xTime, { highPass, lowPass, bandPass }); 
		yTime = xTime - bandPass;
	}
	void ProcessOff(Input xTime, Output yTime) { yTime = xTime; }
};

class IIR2ndTimeVaryingBellFilter : public AsynchronousBase<IIR2ndTimeVaryingBellFilter>
{
	friend Base<IIR2ndTimeVaryingBellFilter>;
public:
	StateVariableFilter Filter;
private:
	struct Coefficients
	{
		float SampleRate = 16e3f;
		int NChannelsIn = 2;
		int BufferSize = 128;
		AsynchronousBufferType AsynchronousBuffer = VARIABLE_SIZE;
	} C;

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

	auto InitializeMembers() 
	{
		auto c = Filter.GetCoefficients();
		c.NChannels = C.NChannelsIn; 
		c.BufferSize = C.BufferSize; 
		c.SampleRate = C.SampleRate; 
		return Filter.Initialize(c); 
	}

	void ProcessOn(Input xTime, Output yTime)
	{
		Eigen::ArrayXXf highPass(xTime.rows(), xTime.cols()), lowPass(xTime.rows(), xTime.cols()), bandPass(xTime.rows(), xTime.cols());
		Filter.Process(xTime, { highPass, lowPass, bandPass });
		yTime = xTime + D.Gain * bandPass;
	}
	void ProcessOff(Input xTime, Output yTime) { yTime = xTime; }
};

class IIR2ndTimeVaryingHighShelfFilter : public AsynchronousBase<IIR2ndTimeVaryingHighShelfFilter>
{
	friend Base<IIR2ndTimeVaryingHighShelfFilter>;
public:
	StateVariableFilter Filter;
private:
	struct Coefficients
	{
		float SampleRate = 16e3f;
		int NChannelsIn = 2;
		int BufferSize = 128;
		AsynchronousBufferType AsynchronousBuffer = VARIABLE_SIZE;
	} C;

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

	auto InitializeMembers() 
	{
		auto c = Filter.GetCoefficients(); 
		c.NChannels = C.NChannelsIn;
		c.BufferSize = C.BufferSize;
		c.SampleRate = C.SampleRate; 
		return Filter.Initialize(c); 
	}

	void ProcessOn(Input xTime, Output yTime)
	{
		Eigen::ArrayXXf highPass(xTime.rows(), xTime.cols()), lowPass(xTime.rows(), xTime.cols()), bandPass(xTime.rows(), xTime.cols());
		Filter.Process(xTime, { highPass, lowPass, bandPass });
		yTime = D.Gain * (bandPass + highPass) + lowPass;
	}
	void ProcessOff(Input xTime, Output yTime) { yTime = xTime; }
};

class IIR2ndTimeVaryingLowShelfFilter : public AsynchronousBase<IIR2ndTimeVaryingLowShelfFilter>
{
	friend Base<IIR2ndTimeVaryingLowShelfFilter>;
public:
	StateVariableFilter Filter;
private:
	struct Coefficients
	{
		float SampleRate = 16e3f;
		int NChannelsIn = 2;
		int BufferSize = 128;
		AsynchronousBufferType AsynchronousBuffer = VARIABLE_SIZE;
	} C;

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

	auto InitializeMembers() 
	{
		auto c = Filter.GetCoefficients(); 
		c.NChannels = C.NChannelsIn;
		c.BufferSize = C.BufferSize; 
		c.SampleRate = C.SampleRate;
		return Filter.Initialize(c); 
	}

	void ProcessOn(Input xTime, Output yTime)
	{
		Eigen::ArrayXXf highPass(xTime.rows(), xTime.cols()), lowPass(xTime.rows(), xTime.cols()), bandPass(xTime.rows(), xTime.cols());
		Filter.Process(xTime, { highPass, lowPass, bandPass });
		yTime = D.Gain * (bandPass + lowPass) + highPass;
	}
	void ProcessOff(Input xTime, Output yTime) { yTime = xTime; }
};