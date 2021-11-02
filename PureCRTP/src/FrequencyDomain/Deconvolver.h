#pragma once
#include "../BaseClasses/PureCRTP.h"
#include "../Utilities/HelperFunctions.h"


struct I::Deconvolver
{
	Complex2D xFreq;
	Boolean speechActivity;
};

class Deconvolver : public BaseFrequencyDomain<Deconvolver, I::Deconvolver>
{
	friend BaseFrequencyDomain<Deconvolver, I::Deconvolver>;

	struct Coefficients
	{
		int NBands = 15;
		int NChannels = 4; // number of channels
		int NLow = 3; // start frame for filter
		int N = 20; // Filter length for one channel ( total filter length is (N + NLow) * NChannels)
		int NSkip = 5; // number of frames to skip between updating the filter in one frequency band (saves computations)
		float Reg = 1e-10f; // regularization parameter for linear systems of equations
	} C;

	struct Parameters
	{
		// the time constant can be calculated from beta and NBeta as: -NBeta/(log(Beta)*filterbankRate). For instance, -400 / (log(0.5) * 16000 / 128) = 4.6 seconds
		float Beta = 0.5f; // forgetting factor
		int NBeta = 400; // number of frames between applying forgetting factor		
	} P;

	struct Data
	{
		Eigen::ConjugateGradient<Eigen::MatrixXcf, Eigen::Lower | Eigen::Upper> SolverCG;
		//std::vector<Eigen::MatrixXcf> BufferInput; // Dimensions: C.nChannels X C.N + C.NLow X C.NBands
		Eigen::MatrixXcf BufferInput;			   // Dimensions: C.nChannels *X C.N + C.NLow X C.NBands
		std::vector<Eigen::MatrixXcf> Filter;      // Dimensions: C.NBands    X C.N			X C.NChannels
		std::vector<Eigen::MatrixXcf> RCov;        // Dimensions: C.NBands    X C.N			X C.N
		std::vector<Eigen::MatrixXcf> RCross;      // Dimensions: C.NBands    X C.N			X C.NChannels
		Eigen::VectorXf Regularization;
		int BufferIndex;
		int BufferLength;
		int FilterCount;
		int FilterLength;
		int DownsamplingIndex;
		int ForgetCount;

		void Reset()
		{
			//for (auto &bufferInput : BufferInput) { bufferInput.setZero(); }
			BufferInput.setZero();
			for (auto &filter : Filter) { filter.setZero(); }
			for (auto &rCov : RCov) 
			{ 
				rCov.setZero(); 
				Eigen::VectorXf drCov(rCov.cols());
				drCov.setConstant(0.01f);
				rCov.diagonal() += drCov;
			}
			for (auto &rCross : RCross)
			{
				rCross.setZero();
				Eigen::VectorXf drCross(rCross.rows());
				drCross.setConstant(1e-12f);
				rCross.col(0) += drCross;
			}
			BufferIndex = 0;
			ForgetCount = 0;
			FilterCount = 0;
			DownsamplingIndex = 0;
		}

		bool InitializeMemory(const Coefficients& c)
		{

			FilterLength = c.N * c.NChannels;
			BufferLength = FilterLength + c.NLow*c.NChannels;

			SolverCG.setMaxIterations(c.NChannels);

			Regularization.resize(FilterLength);
			Regularization.setConstant(c.Reg);
			//BufferInput.resize(c.NChannels);
			//for (auto &bufferInput : BufferInput) {
			//	bufferInput.resize(BufferLength, c.NBands);
			//}
			BufferInput.resize(BufferLength, c.NBands);
			
			Filter.resize(c.NBands);
			RCov.resize(c.NBands);
			RCross.resize(c.NBands);
			for (auto i = 0; i < c.NBands; i++){
				Filter[i].resize(FilterLength, c.NChannels);
				RCov[i].resize(FilterLength, FilterLength);
				RCross[i].resize(FilterLength, c.NChannels);
			}

			if (c.NSkip > c.NBands) { return false; } // if we skip more bands than there are then some bands will never be updated

			return true;
		}

		size_t GetAllocatedMemorySize() const
		{
			size_t size = 0;
			//for (auto &bufferInput : BufferInput) {
			//	size += bufferInput.GetAllocatedMemorySize();
			//}
			size += BufferInput.GetAllocatedMemorySize();
			for (auto i = 0; i < Filter.size(); i++) {
				size += Filter[i].GetAllocatedMemorySize();
				size += RCov[i].GetAllocatedMemorySize();
				size += RCross[i].GetAllocatedMemorySize();
			}
			size += Regularization.GetAllocatedMemorySize();
			return size;
		}

		void OnParameterChange(const Parameters& p, const Coefficients& c) { }

	} D;

	void ProcessOn(Input input, Output yFreq)
	{
		//// copy input into D.BufferInput
		//int channelSelect = D.DownsamplingIndex;
		//for (auto ichan = 0; ichan < C.NChannels; ichan++, channelSelect++)
		//{
		//	if (channelSelect >= C.NChannels) { channelSelect = 0; }
		//	D.BufferInput[ichan].row(D.BufferIndex) = input.xFreq.col(channelSelect).transpose();
		//}
		auto len1 = std::min(D.BufferIndex + C.NChannels, D.BufferLength) - D.BufferIndex;
		auto len2 = C.NChannels - len1;
		D.BufferInput.middleRows(D.BufferIndex, len1) = input.xFreq.leftCols(len1).transpose();
		D.BufferInput.topRows(len2) = input.xFreq.rightCols(len2).transpose();

		// copy input to output
		yFreq = input.xFreq;

		int updateFilterBin = D.FilterCount;
		for (auto ibin = 0; ibin < C.NBands; ibin++)
		{
			// current input frame across channels
			const Eigen::VectorXcf xFreq = input.xFreq.row(ibin).transpose().matrix();

			// saved buffer across time and channels
			int indexStart = D.BufferIndex + C.NLow*C.NChannels; 
			if (indexStart >= D.BufferLength) { indexStart -= D.BufferLength; }
			const auto len1 = std::min(D.BufferLength - indexStart, D.FilterLength);
			const auto len2 = D.FilterLength - len1;
			Eigen::VectorXcf BufferInput(D.FilterLength);
			BufferInput.head(len1) = D.BufferInput.col(ibin).segment(indexStart,len1);
			BufferInput.tail(len2) = D.BufferInput.col(ibin).head(len2);

			//assignCircularEigen(BufferInput.data(), &D.BufferInput(0,ibin), indexStart, len1, len2);

			// filter and subtract from output
			yFreq.row(ibin).matrix() -= BufferInput.transpose() * D.Filter[ibin].conjugate();

			if (input.speechActivity)
			{
				const float invTheta = 1.f / std::max(input.xFreq.row(ibin).abs2().mean(), 1e-10f);
				D.RCov[ibin] += invTheta * BufferInput * BufferInput.adjoint();
				D.RCross[ibin] += invTheta * BufferInput * xFreq.adjoint();
			}

			if (D.ForgetCount == P.NBeta)
			{
				D.RCov[ibin] *= P.Beta;
				D.RCross[ibin] *= P.Beta;
			}

			if (updateFilterBin == 0)
			{
				Eigen::MatrixXcf cov = D.RCov[ibin];
				cov += D.Regularization.asDiagonal();
				// solve linear system of equations: cov * D.Filter[ibin] = D.RCross[ibin]
				D.Filter[ibin] = D.SolverCG.compute(cov).solveWithGuess(D.RCross[ibin], D.Filter[ibin]);
			}
			updateFilterBin++;
			if (updateFilterBin >= C.NSkip) { updateFilterBin = 0; }
		}

		yFreq = (yFreq.abs2() < input.xFreq.abs2()).select(yFreq, input.xFreq);

		D.DownsamplingIndex--;
		if (D.DownsamplingIndex < 0) { D.DownsamplingIndex = C.NChannels - 1; }
		
		D.BufferIndex = D.BufferIndex - C.NChannels;
		if (D.BufferIndex < 0) { D.BufferIndex = D.BufferLength - 1; }

		D.ForgetCount++;
		if (D.ForgetCount == P.NBeta) { D.ForgetCount = 0; }

		D.FilterCount++;
		if (D.FilterCount == C.NSkip) { D.FilterCount = 0; }

	}

	void ProcessOff(Input input, Output yFreq) { yFreq = input.xFreq; }
};