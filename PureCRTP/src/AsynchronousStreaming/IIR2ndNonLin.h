#pragma once
#include "../BaseClasses/PureCRTP.h"
#include "../Utilities/ApproximationMath.h"

class IIR2ndNonLinLowpass : public Base<IIR2ndNonLinLowpass>
{
	friend Base<IIR2ndNonLinLowpass>;

	struct Coefficients
	{
		int NChannels = 2;
		float SampleRate = 16e3f;
	} C;

	struct Parameters
	{
		ConstrainedType<float> Resonance = { static_cast<float>(BUTTERWORTH_Q), 1e-3f, 100.f };
		ConstrainedType<float> Cutoff = { 1000.f, 1.f, 20000.f }; // the max Cutoff should really depend on the samplerate
		ConstrainedType<int> Iterations = { 4, 1, 50 };
	} P;

	struct Data
	{
		Eigen::ArrayXf XOld, Y1Z, Y2Z, Tanh2, Tanh1;
		float g;
		void Reset()
		{
			XOld.setZero();
			Y1Z.setZero();
			Y2Z.setZero();
			Tanh2.setZero();
			Tanh1.setZero();
		}
		bool InitializeMemory(const Coefficients& c)
		{
			XOld.resize(c.NChannels);
			Y1Z.resize(c.NChannels);
			Y2Z.resize(c.NChannels);
			Tanh2.resize(c.NChannels);
			Tanh1.resize(c.NChannels);
			return true;
		}
		size_t GetAllocatedMemorySize() const
		{
			auto size = XOld.GetAllocatedMemorySize();
			size += Y1Z.GetAllocatedMemorySize();
			size += Y2Z.GetAllocatedMemorySize();
			size += Tanh2.GetAllocatedMemorySize();
			size += Tanh1.GetAllocatedMemorySize();
			return size;

		}
		void OnParameterChange(const Parameters& p, const Coefficients& c)
		{
			g = std::min(0.75f, TanApprox(static_cast<float>(PI) * p.Cutoff / (2 * c.SampleRate))); // hardcoded to 2x oversampling so samplerate is 2*c.SampleRate. Limited to avoid numerical problems.
		}
	} D;

	void ProcessOn(Input xTime, Output yTime)
	{
		for (auto sample = 0; sample < xTime.rows(); sample++)
		{
			for (auto channel = 0; channel < xTime.cols(); channel++) // for some strange reason, profiling is faster when iterating accross column before row...
			{
				float v1Iter, v2Iter;
				// hardcoded to 2x oversampling with simple average as anti-aliasing filter
				// ---- average of last and current sample ----
				auto xIn = .5f*(xTime(sample, channel) + D.XOld(channel));
				D.XOld(channel) = xTime(sample, channel);
				xIn = TanhApprox(xIn);
				for (auto j = 0; j < P.Iterations; j++)
				{
					v2Iter = D.g * (D.Tanh1(channel) - D.Tanh2(channel)) + D.Y2Z(channel);
					D.Tanh2(channel) = TanhApprox(v2Iter);
					v1Iter = D.g * (xIn - D.Tanh1(channel)) + P.Resonance * D.Tanh2(channel) + D.Y1Z(channel);
					D.Tanh1(channel) = TanhApprox(v1Iter);
				}
				D.Y1Z(channel) = 2.f * (v1Iter - P.Resonance * v2Iter) - D.Y1Z(channel);
				D.Y2Z(channel) = 2.f * v2Iter - D.Y2Z(channel);
				const auto yOld = v2Iter; // save output for downsampling

				// ----- process current sample ----
				xIn = TanhApprox(xTime(sample, channel));
				for (auto j = 0; j < P.Iterations; j++)
				{
					v2Iter = D.g * (D.Tanh1(channel) - D.Tanh2(channel)) + D.Y2Z(channel);
					D.Tanh2(channel) = TanhApprox(v2Iter);
					v1Iter = D.g * (xIn - D.Tanh1(channel)) + P.Resonance * D.Tanh2(channel) + D.Y1Z(channel);
					D.Tanh1(channel) = TanhApprox(v1Iter);
				}
				D.Y1Z(channel) = 2.f * (v1Iter - P.Resonance * v2Iter) - D.Y1Z(channel);
				D.Y2Z(channel) = 2.f * v2Iter - D.Y2Z(channel);
				// downsample to original samplerate using simple average
				yTime(sample, channel) = .5f*(v2Iter + yOld);
			}
		}
	}

	void ProcessOff(Input xTime, Output yTime) { yTime = xTime; }
};

class IIR2ndNonLinBandpass : public Base<IIR2ndNonLinBandpass>
{
	friend Base<IIR2ndNonLinBandpass>;

	struct Coefficients
	{
		int NChannels = 2;
		float SampleRate = 16e3f;
	} C;

	struct Parameters
	{
		ConstrainedType<float> Resonance = { static_cast<float>(BUTTERWORTH_Q), 1e-3f, 100.f };
		ConstrainedType<float> Cutoff = { 1000.f, 1.f, 20000.f }; // the max Cutoff should really depend on the samplerate
		ConstrainedType<int> Iterations = { 4, 1, 50 };
	} P;

	struct Data
	{
		Eigen::ArrayXf XOld, Y1Z, Y2Z, Y3Z, Tanh3, Tanh2, Tanh1;
		float g;
		void Reset()
		{
			XOld.setZero();
			Y1Z.setZero();
			Y2Z.setZero();
			Y3Z.setZero();
			Tanh3.setZero();
			Tanh2.setZero();
			Tanh1.setZero();
		}
		bool InitializeMemory(const Coefficients& c)
		{
			XOld.resize(c.NChannels);
			Y1Z.resize(c.NChannels);
			Y2Z.resize(c.NChannels);
			Y3Z.resize(c.NChannels);
			Tanh3.resize(c.NChannels);
			Tanh2.resize(c.NChannels);
			Tanh1.resize(c.NChannels);
			return true;
		}
		size_t GetAllocatedMemorySize() const
		{
			auto size = XOld.GetAllocatedMemorySize();
			size += Y1Z.GetAllocatedMemorySize();
			size += Y2Z.GetAllocatedMemorySize();
			size += Y3Z.GetAllocatedMemorySize();
			size += Tanh3.GetAllocatedMemorySize();
			size += Tanh2.GetAllocatedMemorySize();
			size += Tanh1.GetAllocatedMemorySize();
			return size;

		}
		void OnParameterChange(const Parameters& p, const Coefficients& c)
		{
			g = std::min(0.75f, TanApprox(static_cast<float>(PI) * p.Cutoff / (2 * c.SampleRate))); // hardcoded to 2x oversampling so samplerate is 2*c.SampleRate. Limited to avoid numerical problems.
		}
	} D;

	void ProcessOn(Input xTime, Output yTime)
	{
		for (auto sample = 0; sample < xTime.rows(); sample++)
		{
			for (auto channel = 0; channel < xTime.cols(); channel++) // for some strange reason, profiling is faster when iterating accross column before row...
			{
				float v1Iter = .0f, v2Iter = .0f, v3Iter = .0f;
				// hardcoded to 2x oversampling with simple average as anti-aliasing filter
				// ---- average of last and current sample ----
				auto xIn = .5f*(xTime(sample, channel) + D.XOld(channel));
				D.XOld(channel) = xTime(sample, channel);
				xIn = TanhApprox(xIn);
				for (auto j = 0; j < P.Iterations; j++)
				{
					v2Iter = D.g * (D.Tanh1(channel) - D.Tanh2(channel)) + D.Y2Z(channel);
					D.Tanh2(channel) = TanhApprox(v2Iter);
					v1Iter = D.g * (xIn - D.Tanh1(channel)) + P.Resonance * D.Tanh2(channel) + D.Y1Z(channel);
					D.Tanh1(channel) = TanhApprox(v1Iter);
					v3Iter = D.Tanh1(channel) - D.g * D.Tanh3(channel) - D.Y3Z(channel);
					D.Tanh3(channel) = TanhApprox(v3Iter);

				}
				D.Y1Z(channel) = 2.f * (v1Iter - P.Resonance * v2Iter) - D.Y1Z(channel);
				D.Y2Z(channel) = 2.f * v2Iter - D.Y2Z(channel);
				D.Y3Z(channel) += 2.f * D.g * v3Iter;
				const auto yOld = v3Iter; // save output for downsampling

				// ----- process current sample ----
				xIn = TanhApprox(xTime(sample, channel));
				for (auto j = 0; j < P.Iterations; j++)
				{
					v2Iter = D.g * (D.Tanh1(channel) - D.Tanh2(channel)) + D.Y2Z(channel);
					D.Tanh2(channel) = TanhApprox(v2Iter);
					v1Iter = D.g * (xIn - D.Tanh1(channel)) + P.Resonance * D.Tanh2(channel) + D.Y1Z(channel);
					D.Tanh1(channel) = TanhApprox(v1Iter);
					v3Iter = D.Tanh1(channel) - D.g * D.Tanh3(channel) - D.Y3Z(channel);
					D.Tanh3(channel) = TanhApprox(v3Iter);
				}
				D.Y1Z(channel) = 2.f * (v1Iter - P.Resonance * v2Iter) - D.Y1Z(channel);
				D.Y2Z(channel) = 2.f * v2Iter - D.Y2Z(channel);
				D.Y3Z(channel) += 2.f * D.g * v3Iter;
				// downsample to original samplerate using simple average
				yTime(sample, channel) = .5f*(v3Iter + yOld);
			}
		}
	}

	void ProcessOff(Input xTime, Output yTime) { yTime = xTime; }
};

class IIR2ndNonLinLowpassStreaming : public AsynchronousStreaming<IIR2ndNonLinLowpassStreaming, IIR2ndNonLinLowpass>
{
	friend AsynchronousStreaming<IIR2ndNonLinLowpassStreaming, IIR2ndNonLinLowpass>;

	int CalculateLatencySamples() const { return 0; }
	int CalculateNChannelsOut(const int nChannels) const { return nChannels; }

	int UpdateCoefficients(decltype(Algo.GetCoefficients())& c, const float sampleRate, const int nChannels, const int bufferSizeExpected, std::vector<int> bufferSizesSuggested) const
	{
		c.NChannels = nChannels;
		c.SampleRate = sampleRate;
		return bufferSizeExpected;
	}
};

class IIR2ndNonLinBandpassStreaming : public AsynchronousStreaming<IIR2ndNonLinBandpassStreaming, IIR2ndNonLinBandpass>
{
	friend AsynchronousStreaming<IIR2ndNonLinBandpassStreaming, IIR2ndNonLinBandpass>;

	int CalculateLatencySamples() const { return 0; }
	int CalculateNChannelsOut(const int nChannels) const { return nChannels; }

	int UpdateCoefficients(decltype(Algo.GetCoefficients())& c, const float sampleRate, const int nChannels, const int bufferSizeExpected, std::vector<int> bufferSizesSuggested) const
	{
		c.NChannels = nChannels;
		c.SampleRate = sampleRate;
		return bufferSizeExpected;
	}
};