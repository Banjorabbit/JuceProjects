#pragma once
#include "../BaseClasses/PureCRTP.h"

class Distortion5thOrderOdd : public AsynchronousBase<Distortion5thOrderOdd>
{
	friend Base<Distortion5thOrderOdd>;

	struct Coefficients 
	{
		int NChannelsIn = 1;
		int BufferSize = 128;
		AsynchronousBufferType AsynchronousBuffer = VARIABLE_SIZE;
	} C;

	struct Parameters 
	{
		float gain = 1.f;
	} P;

	struct Data
	{
		void Reset() {	}
		bool InitializeMemory(const Coefficients& c) { return true; }
		size_t GetAllocatedMemorySize() const { return 0; }
		void OnParameterChange(const Parameters& p, const Coefficients& c) { }
	} D;

	void ProcessOn(Input xTime, Output yTime) 
	{
		const Eigen::ArrayXXf y = (xTime.block(0,0, C.BufferSize, C.NChannelsIn) * P.gain).min(1.f).max(-1.f) * 1.875f;
		const Eigen::ArrayXXf y2 = y * y;
		const Eigen::ArrayXXf y3 = y * y2;
		yTime = y + (0.016181728395062f * y2 - 0.189629629629630f) * y3;
	}

	void ProcessOff(Input xTime, Output yTime) { yTime.setZero(); }
};
