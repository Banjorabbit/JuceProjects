#pragma once
#include "../BaseClasses/PureCRTP.h"

class SeparateTransient : public AsynchronousBase<SeparateTransient>
{
	friend Base<SeparateTransient>;

	struct Coefficients 
	{
		float SampleRate = 44100.f;
		int BufferSize = 128;
		int NChannelsIn = 2;
		AsynchronousBufferType AsynchronousBuffer = VARIABLE_SIZE;
	} C;
	struct Parameters 
	{
		float AbsTConstant = 0.1f;
		float Threshold = 4.f;
	} P;

	struct Data 
	{
		Eigen::ArrayXf XShort, XLong;
		float AbsLambda;
		void Reset() 
		{
			XShort.setZero();
			XLong.setZero();
		}
		bool InitializeMemory(const Coefficients& c)
		{
			XShort.resize(c.NChannelsIn);
			XLong.resize(c.NChannelsIn);
			return true;
		}
		size_t GetAllocatedMemorySize() const 
		{
			size_t size = XShort.GetAllocatedMemorySize();
			size += XLong.GetAllocatedMemorySize();
			return size; 
		}
		void OnParameterChange(const Parameters& p, const Coefficients& c) 
		{
			AbsLambda = 1.f - expf(-1.f / (c.SampleRate * p.AbsTConstant));
		}
	} D;

	void ProcessOn(Input xTime, Output yTime) 
	{
		for (auto channel = 0; channel < C.NChannelsIn; channel++)
		{
			for (auto i = 0; i < C.BufferSize; i++)
			{
				float xAbs = std::fabs(xTime(i, channel));
				D.XLong(channel) += D.AbsLambda * (xAbs - D.XLong(channel));
				float xMax = std::max(xAbs, D.XShort(channel));
				D.XShort(channel) = xMax + D.AbsLambda * (xAbs - xMax);
				yTime(i, channel) = std::min(1.f, std::max(0.f, D.XShort(channel) - D.XLong(channel) * P.Threshold) / D.XShort(channel)) * xTime(i, channel);
			}
		}
	}

	void ProcessOff(Input xTime, Output yTime) { yTime = xTime; }
};