#pragma once
#include "BaseClasses/PureCRTP.h"

struct I::CubicSpline
{
	Real2D Input;
	Real2D Output;
	Real2D XDesired;
};

// Implementation of Cardinal spline: https://en.wikipedia.org/wiki/Cubic_Hermite_spline
// Default Parameters (sfactor = 0.f) gives the Catmull-Rom spline.
// Basis functions are hardcoded for efficiency.
//
// author: Kristian Timm Andersen
class CubicSpline : public Base<CubicSpline, I::CubicSpline>
{
	friend Base<CubicSpline, I::CubicSpline>;

	struct Coefficients { } C;

	struct Parameters 
	{
		float sfactor = 0.f;
	} P;

	struct Data {
		void Reset() {}
		bool InitializeMemory(const Coefficients& c) { return true; }
		size_t GetAllocatedMemorySize() const { return 0; }
		void OnParameterChange(const Parameters& p, const Coefficients& c) {}
	} D;

	// class used to sort input frequencies and gain
	class sort_indices
	{
	private:
		const float* mparr;
	public:
		sort_indices(const float* parr) : mparr(parr) {}
		bool operator()(int i, int j) const { return mparr[i] < mparr[j]; }
	};

	void ProcessOn(Input x, Output y) 
	{
		// for each spline
		for (auto channel = 0; channel < x.Input.cols(); channel++)
		{
			// sort in ascending order according to x.Input
			Eigen::ArrayXi indexSorted = Eigen::ArrayXi::LinSpaced(static_cast<int>(x.Input.rows()), 0, static_cast<int>(x.Input.rows()) - 1);
			std::sort(indexSorted.data(), indexSorted.data() + indexSorted.size(), sort_indices(x.Input.data()));
			Eigen::ArrayXf input = x.Input(indexSorted, channel);
			Eigen::ArrayXf output = x.Output(indexSorted, channel);

			// mk is tangent in starting points
			Eigen::ArrayXf mk(input.rows());
			auto d = input.rows() - 2;
			mk.segment(1, d) = (1 - P.sfactor) * (output.segment(2,d) - output.segment(0,d)) / (input.segment(2,d) - input.segment(0,d));
			mk.head(1) = (output(1) - output(0)) / (input(1) - input(0));
			mk.tail(1) = (output(output.rows()-1) - output(output.rows()-2)) / (input(input.rows()-1) - input(input.rows()-2));

			// find index of interpolated points where they cross starting points
			int I = 1;
			Eigen::ArrayXi index(input.rows());
			index(0) = 0;
			for (int i = 0; i < x.XDesired.rows(); i++)
			{
				if (x.XDesired(i, channel) > input(I))
				{
					index(I) = i;
					I += 1;
				}
			}
			index.tail(index.rows()-I) = static_cast<int>(x.XDesired.rows());

			// interpolation
			for (auto i = 0; i < input.rows() - 1; i++)
			{
					float H = input(i + 1) - input(i);
					// T is interpolated points between two starting points normalized to the interval (0,1)
					Eigen::ArrayXf T = (x.XDesired.block(index(i), channel, index(i + 1) - index(i), 1) - input(i)) / H;
					// Interpolation function written as a polynomial of T, i.e. y = x + T*(x1 + T*(x2 + T*(x3)));
					y.block(index(i), channel, index(i + 1) - index(i), 1) = output(i) + T * (H*mk(i) + T * ((3.f * (output(i + 1) - output(i)) - H * (2.f * mk(i) + mk(i + 1))) + T * (2.f * (output(i) - output(i + 1)) + H * (mk(i) + mk(i + 1)))));
			}
		}
	}

	void ProcessOff(Input x, Output y) { y.setZero(); }
};