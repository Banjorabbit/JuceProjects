#include "EchoCancellerNLMS.h"
#include "../Utilities/HelperFunctions.h"

using namespace Eigen;
using namespace std;

void EchoCancellerNLMS::ProcessOn(Input xFreq, Output yFreq)
{
	D.CircCounter = D.CircCounter <= 0 ? C.FilterLength - 1 : D.CircCounter - 1;
	// Put loopback signals into buffer
	D.BuffersLoopback.row(D.CircCounter) = xFreq.Loopback.transpose();
		
	// update loopback variance
	D.LoopbackVariance += D.Lambda * (xFreq.Loopback.abs2() - D.LoopbackVariance); // variance of loopback signal i

	yFreq = xFreq.Input;
	Eigen::ArrayXXf pMin = yFreq.abs2();
	//  filter and subtract from input
	const int conv_length1 = C.FilterLength - D.CircCounter;
	for (auto channel = 0; channel < C.NChannels; channel++)
	{
		// estimated transferfunction x loopback
		for (auto ibin = 0; ibin < C.NBands; ibin++)
		{
			std::complex<float> micEst = 0.f;
			for (auto i = 0; i < conv_length1; i++)
			{
				micEst += D.BuffersLoopback(D.CircCounter + i, ibin) * D.Filters[channel](i, ibin);
			}
			for (auto i = 0; i < D.CircCounter; i++)
			{
				micEst += D.BuffersLoopback(i, ibin) * D.Filters[channel](conv_length1 + i, ibin);
			}
			std::complex<float> error = xFreq.Input(ibin, channel) - micEst;
			float pNew = norm(error);
			if (pNew < pMin(ibin, channel)) {
				yFreq(ibin, channel) = error;
				pMin(ibin, channel) = pNew;
			}

			const float pNewLimit = norm(micEst) * D.NearendLimit;
			D.NearendVariance(ibin, channel) += D.Lambda * (std::max(pNew, pNewLimit) - D.NearendVariance(ibin, channel));

			const float p = D.Momentums(ibin, channel) + C.FilterLength * D.CoefficientVariance(ibin, channel);
			float q = p / (C.FilterLength * D.NearendVariance(ibin, channel) + (C.FilterLength + 2) * p * D.LoopbackVariance(ibin) + 1e-30f);
			q = std::min(q, 1.f / (D.LoopbackVariance(ibin) * C.FilterLength + 1e-30f));
			D.Momentums(ibin, channel) = std::max((1.f - q * D.LoopbackVariance(ibin)) * p, 0.01f);
			const std::complex<float> W = q * error;

			const std::complex<float> temp1 = W * std::conj(D.BuffersLoopback(D.CircCounter, ibin));
			D.Filters[channel](0, ibin) += temp1;
			for (auto i = 1; i < conv_length1; i++)
			{
				D.Filters[channel](i, ibin) += W * std::conj(D.BuffersLoopback(D.CircCounter + i, ibin));
			}
			for (auto i = 0; i < D.CircCounter; i++)
			{
				D.Filters[channel](conv_length1 + i, ibin) += W * std::conj(D.BuffersLoopback(i, ibin));
			}
			// update coefficient variance
			D.CoefficientVariance(ibin, channel) += D.Lambda * (temp1.real()*temp1.real() + temp1.imag()*temp1.imag() - D.CoefficientVariance(ibin, channel));
		}
	}
}

void EchoCancellerNLMS::ProcessOff(Input xFreq, Output yFreq) { yFreq = xFreq.Input; }
