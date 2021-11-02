#include "EchoCancellerNLMS.h"
#include "../Utilities/HelperFunctions.h"
// ref implementation 11-18us
// 20-32

constexpr int NBANDS = 257;

using namespace Eigen;
void EchoCancellerNLMS::ProcessOn(Input xFreq, Output yFreq)
{
	D.CircCounter = D.CircCounter <= 0 ? C.FilterLength - 1 : D.CircCounter - 1;
	// Put loopback signals into buffer
	//D.BuffersLoopback.col(D.CircCounter) = xFreq.Loopback;
	auto *p1 = const_cast<std::complex<float>*>(xFreq.Loopback.data());
	assignEigen(&D.BuffersLoopback(0, D.CircCounter), p1, C.NBands);

	// update loopback variance
	//D.LoopbackVariance += D.Lambda * (xFreq.Loopback.abs2() - D.LoopbackVariance); // variance of loopback signal i
	float temp[NBANDS], temp2[NBANDS];
	assignAbs2Eigen(temp, p1, C.NBands);
	subAssignEigen(temp, D.LoopbackVariance.data(), C.NBands);
	multScalarAssignEigen(temp, D.Lambda, C.NBands);
	addAssignEigen(D.LoopbackVariance.data(), temp, C.NBands);

	//yFreq = xFreq.Input;
	p1 = const_cast<std::complex<float>*>(xFreq.Input.data());
	assignEigen(yFreq.data(), p1, C.NBands * C.NChannels);
	
	//  filter and subtract from input
	std::complex<float> ctemp[NBANDS], ctemp2[NBANDS];
	float *pMin = temp;
	std::complex<float> *micEst = ctemp;

	const int conv_length1 = C.FilterLength - D.CircCounter;
	ArrayXcf error(C.NBands), W(C.NBands), temp3(C.NBands);
	ArrayXf micEstPow(C.NBands), errorPow(C.NBands), p(C.NBands), q(C.NBands);
	Array<bool, Dynamic, 1> pMinSelect(C.NBands);
	for (auto channel = 0; channel < C.NChannels; channel++)
	{
		//pMin = yFreq.col(channel).abs2();
		assignAbs2Eigen(pMin, &yFreq(0, channel), C.NBands);

		//micEst = D.BuffersLoopback.col(D.CircCounter) * D.Filters[channel].col(0);
		//for (auto i = 1; i < conv_length1; i++)//, p2 += C.NBands, p3 += C.NBands)
		//{
		//	micEst += D.BuffersLoopback.col(D.CircCounter + i) * D.Filters[channel].col(i);
		//}
		//for (auto i = 0; i < D.CircCounter; i++)//, p2 += C.NBands, p3 += C.NBands)
		//{
		//	micEst += D.BuffersLoopback.col(i) * D.Filters[channel].col(conv_length1 + i);
		//}
		std::complex<float> *cp1 = micEst;
		std::complex<float> *cp2 = &D.BuffersLoopback(0,D.CircCounter);
		std::complex<float> *cp3 = &D.Filters[channel](0, 1);
		assignMultEigen(cp1, cp2, cp3, C.NBands);
		cp2 += C.NBands;
		for (auto i = 1; i < conv_length1; i++, cp2 += C.NBands, cp3 += C.NBands)
		{
			addAssignMultEigen(cp1, cp2, cp3, C.NBands);
		}
		cp2 = D.BuffersLoopback.data();
		for (auto i = 0; i < D.CircCounter; i++, cp2 += C.NBands, cp3 += C.NBands)
		{
			addAssignMultEigen(cp1, cp2, cp3, C.NBands);
		}

		//micEstPow = micEst.abs2();
		assignAbs2Eigen(micEstPow.data(), micEst, C.NBands);

		//error = xFreq.Input.col(channel) - micEst;
		cp1 = const_cast<std::complex<float>*>(&xFreq.Input(0, channel));
		assignSubEigen(error.data(), cp1, micEst, C.NBands);

		//errorPow = error.abs2();
		assignAbs2Eigen(errorPow.data(), error.data(), C.NBands);

		//pMinSelect = errorPow < pMin;
		//yFreq.col(channel) = pMinSelect.select(error, yFreq.col(channel));
		//pMin = pMinSelect.select(errorPow, pMin);
		float *p1 = errorPow.data();
		float *p2 = pMin;
		cp1 = &yFreq(0, channel);
		cp2 = error.data();
		for (auto i = 0; i < C.NBands; i++, p1++, p2++, cp1++, cp2++)
		{
			if (*p1 < *p2)
			{
				*cp1 = *cp2;
				*p2 = *p1;
			}
		}
		
		//D.NearendVariance.col(channel) += D.Lambda * (errorPow.max(micEstPow* D.NearendLimit) - D.NearendVariance.col(channel));
		assignMultScalarEigen(temp, micEstPow.data(), D.NearendLimit, C.NBands);
		p1 = temp;
		p2 = errorPow.data();
		for (auto i = 0; i < C.NBands; i++, p1++, p2++) { if (p2 > p1) { p1 = p2; } }
		subAssignEigen(temp, &D.NearendVariance(0,channel), C.NBands);
		multScalarAssignEigen(temp, D.Lambda, C.NBands);
		addAssignEigen(&D.NearendVariance(0, channel), temp, C.NBands);

		//p = D.Momentums.col(channel) + C.FilterLength * D.CoefficientVariance.col(channel);
		p1 = p.data();
		assignMultScalarEigen(p1, &D.CoefficientVariance(0, channel), C.FilterLength, C.NBands);
		addAssignEigen(p1, &D.Momentums(0, channel), C.NBands);

		//q = p / (C.FilterLength * D.NearendVariance.col(channel) + (C.FilterLength + 2) * p * D.LoopbackVariance + 1e-30f);
		p1 = q.data();
		assignMultScalarEigen(p1, p.data(), (C.FilterLength + 2), C.NBands);
		multAssignEigen(p1, D.LoopbackVariance.data(), C.NBands);
		addScalarAssignEigen(p1, 1e-30f, C.NBands);
		assignMultScalarEigen(temp, &D.NearendVariance(0, channel), C.FilterLength, C.NBands);
		addAssignEigen(temp, p1, C.NBands);
		assignDivEigen(p1, p.data(), temp, C.NBands);

		//q = q.min( 1.f / (D.LoopbackVariance * C.FilterLength + 1e-30f));
		assignMultScalarEigen(temp, D.LoopbackVariance.data(), C.FilterLength, C.NBands);
		addScalarAssignEigen(temp, 1e-30f, C.NBands);
		assignInvEigen(temp2, temp, C.NBands);
		p1 = q.data();
		p2 = temp2;
		for (auto i = 0; i < C.NBands; i++, p1++, p2++)
		{
			if (*p2 < *p1)
			{
				*p1 = *p2;
			}
		}

		//D.Momentums.col(channel) = ((1.f - q * D.LoopbackVariance) * p).max(0.01f);
		assignMultEigen(temp, D.LoopbackVariance.data(), q.data(), C.NBands);
		D.Momentums.col(channel) = 1.f - Map<ArrayXf>(temp,NBANDS);
		multAssignEigen(&D.Momentums(0,channel), p.data(), C.NBands);
		p1 = &D.Momentums(0, channel);
		for (auto i = 0; i < C.NBands; i++, p1++)
		{
			if (*p1 < 0.01f)
			{
				*p1 = 0.01f;
			}
		}

		//W = q * error;
		assignMultEigen(W.data(), q.data(), error.data(), C.NBands);

		//temp = W * D.BuffersLoopback.col(D.CircCounter).conjugate();
		assignMultConjEigen(temp3.data(), W.data(), &D.BuffersLoopback(0, D.CircCounter), C.NChannels);

		//D.Filters[channel].col(0) += temp;
		//for (auto i = 1; i < conv_length1; i++)
		//{
		//	D.Filters[channel].col(i) += W * D.BuffersLoopback.col(D.CircCounter + i).conjugate();
		//}
		//for (auto i = 0; i < D.CircCounter; i++)
		//{
		//	D.Filters[channel].col(conv_length1 + i) += W * D.BuffersLoopback.col(i).conjugate();
		//}
		cp1 = D.Filters[channel].data();
		addAssignEigen(cp1, temp3.data(), C.NBands);
		cp1 += C.NBands;
		cp2 = &D.BuffersLoopback(0, D.CircCounter + 1);
		for (auto i = 1; i < conv_length1; i++, cp1 += C.NBands, cp2 += C.NBands)
		{
			addAssignMultConjEigen(cp1, W.data(), cp2, C.NBands);
		}
		cp2 = D.BuffersLoopback.data();
		for (auto i = 0; i < D.CircCounter; i++, cp1 += C.NBands, cp2 += C.NBands)
		{
			addAssignMultConjEigen(cp1, W.data(), cp2, C.NBands);
		}

		// update coefficient variance
		//D.CoefficientVariance.col(channel) += D.Lambda * (temp3.abs2() - D.CoefficientVariance.col(channel));
		assignAbs2Eigen(temp, temp3.data(), C.NBands);
		subAssignEigen(temp, &D.CoefficientVariance(0, channel), C.NBands);
		multScalarAssignEigen(temp, D.Lambda, C.NBands);
		addAssignEigen(&D.CoefficientVariance(0, channel), temp, C.NBands);
	}
}

void EchoCancellerNLMS::ProcessOff(Input xFreq, Output yFreq) { yFreq = xFreq.Input; }