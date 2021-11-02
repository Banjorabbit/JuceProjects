#pragma once
#include "../BaseClasses/PureCRTP.h"
#include "EchoCancellerNLMS.h"

// predictor class that uses the NLMS echo canceller to predict the input.
// The output is the error of the prediction.
class PredictorNLMS : public BaseFrequencyDomain<PredictorNLMS>
{
	friend BaseFrequencyDomain<PredictorNLMS>;

	VectorAlgo<EchoCancellerNLMS> Predictor;

private:

	struct Coefficients
	{
		float FilterbankRate = 125.f;
		int FilterOrder = 4;
		int NBands = 257;
		int NChannels = 2;
	} C;

	struct Parameters
	{
	} P;

	struct Data
	{
		Eigen::ArrayXXcf InputOld;
		void Reset();
		bool InitializeMemory(const Coefficients& c);
		size_t GetAllocatedMemorySize() const;
		void OnParameterChange(const Parameters& p, const Coefficients& c);
	} D;

	DEFINEMEMBERALGORITHMS(1, Predictor);

	bool InitializeMembers();

	void ProcessOn(Input xFreq, Output yFreq);

	void ProcessOff(Input xFreq, Output yFreq);
};