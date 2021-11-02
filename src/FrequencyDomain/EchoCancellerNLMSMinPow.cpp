#include "EchoCancellerNLMSMinPow.h"

using namespace Eigen;
using namespace std;

void EchoCancellerNLMSMinPow::ProcessOn(Input xFreq, Output yFreq)
{

	D.CircCounter = D.CircCounter <= 0 ? (int)D.BuffersLoopback.rows() - 1 : D.CircCounter - 1;
	// update circular loopback buffer
	D.BuffersLoopback.row(D.CircCounter) = xFreq.Loopback.transpose();

	// update loopback variance
	for (auto ibin = 0; ibin < C.NBands; ibin++)
	{
		D.LoopbackVariance.col(ibin) += D.Lambda * (norm(xFreq.Loopback(ibin)) - D.LoopbackVariance.col(ibin));
	}


	// copy input to output
	yFreq = xFreq.Input;
	ArrayXXf powMin = xFreq.Input.abs2();
	// calculate filter lengths
	for (auto ichan = 0; ichan < C.NChannels; ichan++)
	{
		for (auto ibin = 0; ibin < C.NBands; ibin++)
		{
			for (auto ifilt = 0; ifilt < C.NFilters; ifilt++)
			{
				const int convLength1 = min(C.FilterLength[ifilt], (int)D.BuffersLoopback.rows() - D.CircCounter);
				const int convLength2 = C.FilterLength[ifilt] - convLength1;
				complex<float> micEst = 0.f;
				for (auto i = 0; i < convLength1; i++)
				{
					micEst += D.BuffersLoopback(D.CircCounter + i, ibin) * D.Filters[ichan][ifilt](i, ibin);
				}
				for (auto i = 0; i < convLength2; i++)
				{
					micEst += D.BuffersLoopback(i, ibin) * D.Filters[ichan][ifilt](convLength1 + i, ibin);
				}
				complex<float> error = xFreq.Input(ibin, ichan) - micEst;
				float powNew = norm(error);
				if (powNew < powMin(ibin, ichan))
				{
					powMin(ibin, ichan) = powNew;
					yFreq(ibin, ichan) = error;
				}

				const float pNewLimit = norm(micEst) * D.NearendLimit;
				D.NearendVariance[ichan](ifilt, ibin) += D.Lambda(ifilt) * (std::max(powNew, pNewLimit) - D.NearendVariance[ichan](ifilt, ibin));

				const float p = D.Momentums[ichan](ifilt, ibin) + C.FilterLength[ifilt] * D.CoefficientVariance[ichan](ifilt, ibin);
				float q = p / (C.FilterLength[ifilt] * D.NearendVariance[ichan](ifilt, ibin) + (C.FilterLength[ifilt] + 2) * p * D.LoopbackVariance(ifilt, ibin) + 1e-10f);
				q = min(q, 1.f / (D.LoopbackVariance(ifilt, ibin) * C.FilterLength[ifilt] + 1e-20f));
				D.Momentums[ichan](ifilt, ibin) = max((1.f - q * D.LoopbackVariance(ifilt, ibin)) * p, 0.01f);
				const complex<float> W = q * error;

				const std::complex<float> temp1 = W * conj(D.BuffersLoopback(D.CircCounter, ibin));
				D.Filters[ichan][ifilt](0, ibin) += temp1;
				for (auto i = 1; i < convLength1; i++)
				{
					D.Filters[ichan][ifilt](i, ibin) += W * conj(D.BuffersLoopback(D.CircCounter + i, ibin));
				}
				for (auto i = 0; i < convLength2; i++)
				{
					D.Filters[ichan][ifilt](convLength1 + i, ibin) += W * conj(D.BuffersLoopback(i, ibin));
				}

				// update coefficient variance
				D.CoefficientVariance[ichan](ifilt, ibin) += D.Lambda(ifilt) * (norm(temp1) - D.CoefficientVariance[ichan](ifilt, ibin));
			}
		}
	}

}

void EchoCancellerNLMSMinPow::ProcessOff(Input xFreq, Output yFreq) { yFreq = xFreq.Input; }

