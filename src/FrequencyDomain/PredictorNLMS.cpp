#include "PredictorNLMS.h"


using namespace Eigen;
using namespace std;

void PredictorNLMS::ProcessOn(Input xFreq, Output yFreq)
{
	for (auto channel = 0; channel < C.NChannels; channel++)
	{
		Predictor[channel].Process({ xFreq.col(channel), D.InputOld.col(channel) }, yFreq.col(channel));
	}
	D.InputOld = xFreq;
}

void PredictorNLMS::ProcessOff(Input xFreq, Output yFreq)
{
	yFreq = xFreq;
}

void PredictorNLMS::Data::Reset()
{
	InputOld.setZero();
}

bool PredictorNLMS::Data::InitializeMemory(const Coefficients& c)
{
	InputOld.resize(c.NBands, c.NChannels);
	return true;
}

size_t PredictorNLMS::Data::GetAllocatedMemorySize() const
{
	size_t size = 0;
	size += InputOld.GetAllocatedMemorySize();
	return size;
}

void PredictorNLMS::Data::OnParameterChange(const Parameters& p, const Coefficients& c)
{}

bool PredictorNLMS::InitializeMembers()
{
	Predictor.resize(C.NChannels);
	auto flag = true;
	for (auto& pred : Predictor)
	{
		auto cNLMS = pred.GetCoefficients();
		cNLMS.FilterbankRate = C.FilterbankRate;
		cNLMS.FilterLength = C.FilterOrder;
		cNLMS.NBands = C.NBands;
		cNLMS.NChannels = 1;
		flag &= pred.Initialize(cNLMS);
	}
	return flag;
}