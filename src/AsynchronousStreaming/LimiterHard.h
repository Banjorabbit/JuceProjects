#include "../BaseClasses/PureCRTP.h"
#include "../FilterMinMax.h"
#include "../CircBuffer.h"

class LimiterHard : public Base<LimiterHard>
{
	friend Base<LimiterHard>;

public:
	StreamingMax FindMax;
	CircBuffer DelayLine;
	int GetLatencySamples() const { return static_cast<int>(C.LookAheadMS / 1000.f * C.SampleRate); }

private:

	struct Coefficients 
	{
		float SampleRate = 16000.f;
		float LookAheadMS = 8.f;
		float HoldTimeMS = 12.f;
		int BufferSize = 128;
		int NChannels = 2;
	} C;

	struct Parameters 
	{
		float ReleaseTConstant = .2f;
		float PreGain = 2.f;
		float PostGain = 1.f;
	} P;

	struct Data 
	{
		float ReleaseUp, AttackDown;
		int LengthMax, LengthDelay;
		float Gain, GainAttackSmooth;
		void Reset() 
		{
			Gain = 1.f;
			GainAttackSmooth = 1.f;
		}
		bool InitializeMemory(const Coefficients& c) 
		{
			LengthMax = static_cast<int>((c.LookAheadMS + c.HoldTimeMS) / 1000.f * c.SampleRate);
			LengthDelay = static_cast<int>(c.LookAheadMS / 1000.f * c.SampleRate);
			return true;
		}
		size_t GetAllocatedMemorySize() const 
		{
			return 0;
		}
		void OnParameterChange(const Parameters& p, const Coefficients& c)
		{
			ReleaseUp = expf(1.f / (c.SampleRate*p.ReleaseTConstant));
			AttackDown = 1.f - expf(-1.f / (c.SampleRate * c.LookAheadMS / 1000.f / 3.f)); // 1st order filter reaches 95% of final value in 3 time constants
		}
	} D;

	DEFINEMEMBERALGORITHMS(2, FindMax, DelayLine);

	auto InitializeMembers()
	{
		auto c = FindMax.GetCoefficients();
		c.Length = D.LengthMax;
		c.NChannels = 1;
		auto flag =  FindMax.Initialize(c);

		auto cDL = DelayLine.GetCoefficients();
		cDL.DelayLength = D.LengthDelay;
		cDL.NChannels = C.NChannels;
		flag &= DelayLine.Initialize(cDL);

		return flag;
	}

	void ProcessOn(Input xTime, Output yTime) 
	{
		// apply pregain
		Eigen::ArrayXXf xPre = xTime * P.PreGain;

		// delay xPre with LookAhead
		DelayLine.Process(xPre, yTime);

		// find max(abs(xTime)) larger than 1 over LookAhead + HoldTime
		Eigen::ArrayXf xMax(C.BufferSize);
		FindMax.Process(xPre.abs().rowwise().maxCoeff().max(1.f), xMax);

		// calculate desired output gain
		for (auto i = 0; i < C.BufferSize; i++)
		{
			float gain = 1.f / xMax(i);
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

class LimiterHardStreaming : public AsynchronousStreaming<LimiterHardStreaming, LimiterHard>
{
	friend AsynchronousStreaming<LimiterHardStreaming, LimiterHard>;

	int CalculateLatencySamples() const { return Algo.GetLatencySamples(); }
	int CalculateNChannelsOut(const int nChannels) const { return nChannels; }

	int UpdateCoefficients(decltype(Algo.GetCoefficients())& c, const float sampleRate, const int nChannels, const int bufferSizeExpected, std::vector<int> bufferSizesSuggested) const
	{
		c.NChannels = nChannels;
		c.SampleRate = sampleRate;
		c.BufferSize = bufferSizeExpected;
		return c.BufferSize;
	}
};