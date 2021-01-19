#pragma once
#include "BaseClasses/PureCRTP.h"
#include "InterpolationCubic.h"

// 2x upsampling using 4 point, 3rd order Hermite interpolation based on :
// http ://yehar.com/blog/wp-content/uploads/2009/08/deip.pdf
//
// This upsampling assumes that there is an additional point after the last samples that has the value 0.
//
// author: Kristian Timm Andersen

class Upsampling2XCubic : public Base<Upsampling2XCubic>
{
	friend Base<Upsampling2XCubic>;

	struct Coefficients {} C;
	struct Parameters {} P;

	struct Data 
	{
		static constexpr float a0 = -0.062f;
		static constexpr float a1 = 0.5625f;
		void Reset() { }
		bool InitializeMemory(const Coefficients& c) { return true; }
		size_t GetAllocatedMemorySize() const { return 0; }
		void OnParameterChange(const Parameters& p, const Coefficients& c) {}
	} D;

	void ProcessOn(Input xTime, Output yTime) 
	{
		const auto size1 = xTime.rows() - 1;
		const auto size3 = xTime.rows() - 3;
		Eigen::ArrayXXf halfPoint(xTime.rows(), xTime.cols());
		halfPoint.topRows(size1) = D.a1 * (xTime.topRows(size1) + xTime.bottomRows(size1));
		halfPoint.row(size1) = D.a1 * xTime.row(size1);

		halfPoint.row(0) += D.a0 * xTime.row(2);
		halfPoint.middleRows(1, size3) += D.a0 * (xTime.topRows(size3) + xTime.bottomRows(size3));
		halfPoint.bottomRows(2) += D.a0 * xTime.middleRows(size3, 2);
		// use Eigen::Map to alternate between xTime and halfPoint
		Eigen::Map<Eigen::ArrayXXf, 0, Eigen::InnerStride<2>>(&yTime(0, 0), xTime.rows(), xTime.cols()) = xTime;
		Eigen::Map<Eigen::ArrayXXf, 0, Eigen::InnerStride<2>>(&yTime(1, 0), xTime.rows(), xTime.cols()) = halfPoint;
	}

	void ProcessOff(Input xTime, Output yTime) { yTime.setZero(); }
};

// power of 2 upsampling
class UpsamplingPower2Cubic : public Base<UpsamplingPower2Cubic>
{
	friend Base<UpsamplingPower2Cubic>;

public:
	Upsampling2XCubic Upsample2X;

private:
	struct Coefficients {} C;
	struct Parameters 
	{
		ConstrainedType<int> UpsamplingFactorPow2 = { 2, 1, 10 }; // code doesn't support <= 0
	} P;

	struct Data
	{
		void Reset() { }
		bool InitializeMemory(const Coefficients& c) { return true; }
		size_t GetAllocatedMemorySize() const { return 0; }
		void OnParameterChange(const Parameters& p, const Coefficients& c) {}
	} D;
	
	DEFINEMEMBERALGORITHMS(1, Upsample2X);

	auto InitializeMembers() { return Upsample2X.Initialize(); }

	void ProcessOn(Input xTime, Output yTime)
	{
		auto size = xTime.rows();
		Eigen::ArrayXXf data(size * (2 << (P.UpsamplingFactorPow2 - 1)), xTime.cols());
		data.topRows(size) = xTime;
		for (auto i = 0; i < P.UpsamplingFactorPow2;i++)
		{
			Eigen::ArrayXXf output(2 * size, xTime.cols());
			Upsample2X.Process(data.topRows(size), output);
			size *= 2;
			data.topRows(size) = output;
		}
		yTime = data;
	}

	void ProcessOff(Input xTime, Output yTime) { yTime.setZero(); }
};