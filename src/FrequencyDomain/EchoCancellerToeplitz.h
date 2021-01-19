#pragma once
#include "../BaseClasses/PureCRTP.h"
#include "../ToeplitzSolver.h"

struct I::EchoCancellerToeplitz
{
	Complex2D Input;
	Complex Loopback;
};

class EchoCancellerToeplitz : public BaseFrequencyDomain<EchoCancellerToeplitz, I::EchoCancellerToeplitz>
{
	friend BaseFrequencyDomain<EchoCancellerToeplitz, I::EchoCancellerToeplitz>;

	ToeplitzSolver TSolver;

private:
	struct Coefficients
	{
		float FilterbankRate = 125.f;
		float FilterLengthMilliseconds = 250.f;
		int NBands = 257;
		int NChannels = 2;
		int nBuffersInAnalysisFrame = 4;
		int nBuffersInSynthesisFrame = 4;
	} C;

	struct Parameters {} P;

	struct Data
	{
		std::vector<Eigen::MatrixXcf> Filters;
		int FilterLength;
		Eigen::MatrixXcf BuffersLoopback;
		int CircCounter;
		float beta;
		Eigen::MatrixXcf Crx;
		std::vector<Eigen::MatrixXcf> Crxy;
		int FilterUpdatesPerFrame;
		int FilterUpdateCounter;

		void Reset()
		{
			BuffersLoopback.setZero();
			for (auto &filter : Filters) { filter.setZero(); }
			Crx.setZero();
			for (auto &crxy : Crxy) { crxy.setZero(); }
			CircCounter = 0;
			FilterUpdateCounter = 0;
		}
		bool InitializeMemory(const Coefficients& c)
		{
			FilterLength = std::max(static_cast<int>(c.FilterLengthMilliseconds * c.FilterbankRate * 1e-3f) - (c.nBuffersInAnalysisFrame + c.nBuffersInSynthesisFrame - 1), 0) + 1;
			FilterUpdatesPerFrame = static_cast<int>(c.NBands / FilterLength); // update all filters after FilterLength frames
			beta = expf(-1.f / (c.FilterbankRate / FilterLength * c.FilterLengthMilliseconds * 1e-3f)); // time constant of covariance update equal to filter length in seconds
			BuffersLoopback.resize(FilterLength, c.NBands);
			Filters.resize(c.NChannels);
			for (auto &filter : Filters) { filter.resize(FilterLength, c.NBands); }
			Crx.resize(FilterLength, c.NBands);
			Crxy.resize(c.NChannels);
			for (auto &crxy : Crxy) { crxy.resize(FilterLength, c.NBands); }
			return true;
		}
		size_t GetAllocatedMemorySize() const
		{
			size_t size = 0;
			for (auto& filter : Filters) { size += filter.GetAllocatedMemorySize(); }
			size += BuffersLoopback.GetAllocatedMemorySize();
			for (auto& crxy : Crxy) { size += crxy.GetAllocatedMemorySize(); }
			return size;
		}
		void OnParameterChange(const Parameters& p, const Coefficients& c) {}
	} D;

	DEFINEMEMBERALGORITHMS(1, TSolver);

	auto InitializeMembers()
	{
		auto pToeplitz = TSolver.GetParameters();
		pToeplitz.regularization = 1e-4f;
		TSolver.SetParameters(pToeplitz);
		return TSolver.Initialize();
	}

	void ProcessOn(Input xFreq, Output yFreq)
	{
		D.CircCounter = D.CircCounter <= 0 ? D.FilterLength - 1 : D.CircCounter - 1;
		int filterLength1 = D.FilterLength - D.CircCounter;
		// Put loopback signals into buffer
		D.BuffersLoopback.row(D.CircCounter) = xFreq.Loopback.matrix().adjoint();

		yFreq = xFreq.Input;
		// RX covariance
		for (int ibin = 0; ibin < C.NBands; ibin++)
		{
			const auto iLoopback = D.BuffersLoopback.col(ibin);
			D.Crx.col(ibin).head(filterLength1) += xFreq.Loopback(ibin) * iLoopback.segment(D.CircCounter, filterLength1);
			D.Crx.col(ibin).segment(filterLength1, D.CircCounter) += xFreq.Loopback(ibin) * iLoopback.head(D.CircCounter);
			for (auto channel = 0; channel < C.NChannels; channel++)
			{
				// cross RX and TX covariance
				D.Crxy[channel].col(ibin).head(filterLength1) += xFreq.Input(ibin, channel) * iLoopback.segment(D.CircCounter, filterLength1);
				D.Crxy[channel].col(ibin).segment(filterLength1, D.CircCounter) += xFreq.Input(ibin, channel) * iLoopback.head(D.CircCounter);
				// do filter
				const auto iFilter = D.Filters[channel].col(ibin);
				yFreq(ibin, channel) -= iLoopback.segment(D.CircCounter, filterLength1).adjoint() * iFilter.head(filterLength1);
				yFreq(ibin, channel) -= iLoopback.head(D.CircCounter).adjoint() * iFilter.tail(D.CircCounter);
			}	
		}

		for (int i = 0; i < D.FilterUpdatesPerFrame; i++)
		{
			// create right-hand side matrix
			Eigen::MatrixXcf RToep(D.FilterLength, C.NChannels);
			for (auto channel = 0; channel < C.NChannels; channel++)
			{
				RToep.col(channel) = D.Crxy[channel].col(D.FilterUpdateCounter);
			}
			Eigen::ArrayXXcf result(D.FilterLength, C.NChannels);
			TSolver.Process({ D.Crx.col(D.FilterUpdateCounter).conjugate(), RToep }, result);
			for (auto channel = 0; channel < C.NChannels; channel++)
			{
				D.Filters[channel].col(D.FilterUpdateCounter) = result.col(channel);
			}

			// Apply forgetting factor to Crx and Crxy
			D.Crx.col(D.FilterUpdateCounter) *= D.beta;
			for (auto channel = 0; channel < C.NChannels; channel++)
			{
				D.Crxy[channel].col(D.FilterUpdateCounter) *= D.beta;
			}

			D.FilterUpdateCounter++;
			if (D.FilterUpdateCounter >= C.NBands) { D.FilterUpdateCounter = 0; }
		}
	}

	void ProcessOff(Input xFreq, Output yFreq) { yFreq = xFreq.Input; }
};