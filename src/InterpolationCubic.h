#pragma once
#include "BaseClasses/PureCRTP.h"

// 4 point, 3rd order Hermite interpolation based on :
// http ://yehar.com/blog/wp-content/uploads/2009/08/deip.pdf
//
// x.Input assumed to be a size 4 array and x.FractionalDelay is
// in the interval[0, 1[ between x.Input[1] and x.Input[2].
//
// author: Kristian Timm Andersen

struct I::InterpolationCubic
{
	Real4 Input;
	Real FractionalDelay; // FractionalDelay and Output must have same size
};

class InterpolationCubic : public Base<InterpolationCubic, I::InterpolationCubic>
{
	friend Base<InterpolationCubic, I::InterpolationCubic>;

	struct Coefficients {} C;
	struct Parameters {} P;

	struct Data 
	{
		void Reset() { }
		bool InitializeMemory(const Coefficients& c) { return true; }
		size_t GetAllocatedMemorySize() const { return 0; }
		void OnParameterChange(const Parameters& p, const Coefficients& c) {}
	} D;

	void ProcessOn(Input x, Output y) 
	{
		const auto c0 = x.Input(1);
		const auto c1 = .5f*(x.Input(2) - x.Input(0));
		const auto c2 = x.Input(0) - 2.5f*x.Input(1) + 2.f*x.Input(2) - .5f*x.Input(3);
		const auto c3 = .5f*(x.Input(3) - x.Input(0)) + 1.5f*(x.Input(1) - x.Input(2));
		y = ((c3*x.FractionalDelay + c2)*x.FractionalDelay + c1)*x.FractionalDelay + c0;
	}

	void ProcessOff(Input x, Output y) { y.setZero(); }
};

// When the fractional delay is constant in all interpolations, this function is faster.
class InterpolationCubicConstantDelay : public Base<InterpolationCubicConstantDelay, I::Real4X, O::Real>
{
	friend Base<InterpolationCubicConstantDelay, I::Real4X, O::Real>;

	struct Coefficients { } C;

	struct Parameters 
	{
		ConstrainedType<float> FractionalDelay = { 0.5f, 0.f, 1.f };
	} P;
	
	struct Data
	{
		float a0, a1, a2, a3;
		void Reset() { }
		bool InitializeMemory(const Coefficients& c) { return true; }
		size_t GetAllocatedMemorySize() const { return 0; }
		void OnParameterChange(const Parameters& p, const Coefficients& c) 
		{
			const auto d2 = p.FractionalDelay * p.FractionalDelay;
			const auto dm1 = p.FractionalDelay - 1.f;
			const auto dm12 = dm1 *dm1;
			a0 = -dm12 * p.FractionalDelay *.5f;
			a1 = dm12 * (1.f + 2.f * p.FractionalDelay) - d2 * dm1 *.5f;
			a2 = dm12 * p.FractionalDelay *.5f + d2 * (3.f - 2.f * p.FractionalDelay);
			a3 = d2 * dm1 *.5f;
		}
	} D;

	void ProcessOn(Input x, Output y)
	{
		// this has been profiled to be faster than both Eigen matrix multiplications (slowest) and for-loops.
		y = (D.a0 * x.row(0) + D.a1 * x.row(1) + D.a2 * x.row(2) + D.a3 * x.row(3)).transpose();
	}

	void ProcessOff(Input x, Output y) { y.setZero(); }
};