#pragma once
#include "BaseClasses/PureCRTP.h"

struct I::ToeplitzSolver
{
	Complex aToeplitz;
	Complex2D BRightHand;
};

// Toeplitz Solver that can solve multiple right hand sides at once.
//
// Input:
// aToeplitz - 1D complex array that corresponds to the first row in lefthand toeplitz matrix
// BRightHand - 2D complex array where each column is a right hand side vector
//
// author: Kristian Timm Andersen
class ToeplitzSolver : public Base<ToeplitzSolver, I::ToeplitzSolver, O::Complex2D>
{
	friend Base<ToeplitzSolver, I::ToeplitzSolver, O::Complex2D>;

	struct Coefficients {} C;

	struct Parameters 
	{
		float regularization = 1e-4f;
	} P;

	struct Data
	{
		void Reset() { }
		bool InitializeMemory(const Coefficients& c) { return true;	}
		size_t GetAllocatedMemorySize() const {	return 0; }
		void OnParameterChange(const Parameters& p, const Coefficients& c) {}
	} D;

	void ProcessOn(Input x, Output y)
	{
		const auto nFilt = x.BRightHand.rows();
		const auto nChan = x.BRightHand.cols();

		Eigen::ArrayXcf c1(nFilt - 1);
		Eigen::ArrayXcf c2(nFilt - 1);
		std::complex<float> r1 = { 0,0 }, r2 = { 0,0 }, r3 = { 0,0 }, r5 = { 0,0 }, r6 = { 0,0 };

		c1.setZero();
		c2.setZero();
		y.setZero();

		float m_R = x.aToeplitz(0).real();
		r1.real(x.aToeplitz(0).real() + P.regularization * m_R);
		r1.imag(x.aToeplitz(0).imag());

		float recip_r1;
		recip_r1 = 1.0f / (r1.real()*r1.real() + r1.imag()*r1.imag());

		for (int ichan = 0; ichan < nChan; ichan++) {
			y(0, ichan) = x.BRightHand(0, ichan) * std::conj(r1) * recip_r1;
		}

		for (int isub = 1; isub < nFilt; isub++) {
			r5 = std::conj(x.aToeplitz(isub));
			r6 = x.aToeplitz(isub);
			if (isub > 1) {
				c1(isub - 1) = r2;
				for (int i = 0; i <= isub - 2; i++) {
					std::complex<float> R_toep_conj = std::conj(x.aToeplitz(i + 1));
					r5 += R_toep_conj * c1(isub - (i + 1));
					r6 += x.aToeplitz(i + 1) * c2(i);
				}
			}
			recip_r1 = 1.0f / (r1.real() * r1.real() + r1.imag() * r1.imag());
			r2 = -r5 * std::conj(r1) * recip_r1;
			r3 = -r6 * std::conj(r1) * recip_r1;
			r1 += r5 * r3;

			if (isub > 1) {
				r6 = c2(0);

				c2(isub - 1) = 0;

				for (int i = 1; i <= isub - 1; i++) {

					r5 = c2(i);
					c2(i) = c1(i) * r3 + r6;
					c1(i) += r6 * r2;
					r6 = r5;
				}
			}
			c2(0) = r3;

			for (int ichan = 0; ichan < nChan; ichan++) {
				r5 = 0.f;

				for (int i = 0; i <= isub - 1; i++) {
					std::complex<float> R_toep_conj = std::conj(x.aToeplitz(i + 1));
					r5 += R_toep_conj * y(isub - (i + 1), ichan);
				}

				std::complex<float> temp = x.BRightHand(isub, ichan) - r5;

				recip_r1 = 1.0f / (r1.real() * r1.real() + r1.imag() * r1.imag());
				r6 = temp * std::conj(r1) * recip_r1;

				for (int i = 0; i <= isub - 1; i++) {
					y(i, ichan) += c2(i) * r6;
				}
				y(isub, ichan) = r6;
			}
		} 
	}

	void ProcessOff(Input x, Output y) { y.setZero(); }
};
