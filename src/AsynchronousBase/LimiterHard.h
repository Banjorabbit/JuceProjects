#pragma once
#include "../BaseClasses/PureCRTP.h"
#include "../FilterMinMax.h"
#include "../CircBuffer.h"

class LimiterHard : public AsynchronousBase<LimiterHard>
{
	friend Base<LimiterHard>;

public:
	StreamingMax FindMax;
	int GetLatencySamples() const { return static_cast<int>(C.LookAheadMS / 1000.f * C.SampleRate); }

private:

	struct Coefficients 
	{
		float SampleRate = 16000.f;
		float LookAheadMS = 8.f;
		float HoldTimeMS = 12.f;
		int BufferSize = 128;
		int NChannelsIn = 2;
		AsynchronousBufferType AsynchronousBuffer = VARIABLE_SIZE;
	} C;

	struct Parameters 
	{
		float ReleaseTConstant = .2f;
		float PreGain = 2.f;
		float PostGain = 1.f;
		float Threshold = 1.f;
	} P;

	struct Data 
	{
		float ReleaseUp, AttackDown;
		int LengthMax, LengthDelay;
		float Gain, GainAttackSmooth;
		int size1, size2, size3, size4;
		Eigen::ArrayXXf DelayLine;
		void Reset() 
		{
			Gain = 1.f;
			GainAttackSmooth = 1.f;
			DelayLine.setZero();
		}
		bool InitializeMemory(const Coefficients& c) 
		{
			LengthMax = static_cast<int>((c.LookAheadMS + c.HoldTimeMS) / 1000.f * c.SampleRate);
			LengthDelay = static_cast<int>(c.LookAheadMS / 1000.f * c.SampleRate);
			DelayLine.resize(LengthDelay, c.NChannelsIn);
			if (c.BufferSize < LengthDelay)
			{
				size1 = c.BufferSize;
				size2 = 0;
				size3 = LengthDelay - c.BufferSize;
				size4 = c.BufferSize;
			}
			else
			{
				size1 = LengthDelay;
				size2 = c.BufferSize - LengthDelay;
				size3 = 0;
				size4 = LengthDelay;
			}
			AttackDown = 1.f - expf(-1.f / (c.SampleRate * c.LookAheadMS / 1000.f / 3.f)); // 1st order filter reaches 95% of final value in 3 time constants
			return true;
		}
		size_t GetAllocatedMemorySize() const 
		{
			return DelayLine.GetAllocatedMemorySize();
		}
		void OnParameterChange(const Parameters& p, const Coefficients& c)
		{
			ReleaseUp = expf(1.f / (c.SampleRate*p.ReleaseTConstant));
		}
	} D;

	DEFINEMEMBERALGORITHMS(1, FindMax);

	auto InitializeMembers()
	{
		auto c = FindMax.GetCoefficients();
		c.Length = D.LengthMax;
		c.NChannels = 1;
		return FindMax.Initialize(c);
	}

	void ProcessOn(Input xTime, Output yTime) 
	{
		// apply pregain
		Eigen::ArrayXXf xPre = xTime * P.PreGain;

		// delay xPre with LookAhead
		yTime.topRows(D.size1) = D.DelayLine.topRows(D.size1);
		yTime.bottomRows(D.size2) = xPre.topRows(D.size2);
		D.DelayLine.topRows(D.size3) = D.DelayLine.bottomRows(D.size3);
		D.DelayLine.bottomRows(D.size4) = xPre.bottomRows(D.size4);

		// find max(abs(xTime)) larger than Threshold over LookAhead + HoldTime
		Eigen::ArrayXf xMax(C.BufferSize);
		FindMax.Process(xPre.abs().rowwise().maxCoeff().max(P.Threshold), xMax);

		// calculate desired output gain
		for (auto i = 0; i < C.BufferSize; i++)
		{
			float gain = P.Threshold / xMax(i);
			D.GainAttackSmooth += D.AttackDown * (gain*.95f - D.GainAttackSmooth); // 95% of final value is reached after 3 time constants
			D.GainAttackSmooth = std::max(D.GainAttackSmooth, gain); // if no reduction, GainAttackSmooth should be 1
			D.Gain = std::min(D.GainAttackSmooth, D.Gain*D.ReleaseUp);
			yTime.row(i) *= D.Gain * P.PostGain;
		}
	}

	void ProcessOff(Input xTime, Output yTime)
	{
		yTime = xTime;
	}
};