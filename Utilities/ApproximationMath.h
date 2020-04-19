#pragma once
#include <cmath>
#include <algorithm>

const double DOUBLE2_FIX_MAGIC = 6755399441055744.0;   // 0x4338000000000000 = 1.5 *2^52

typedef union {
	double D;
	uint32_t I[2];
} DoubleUnion;

inline int Real2Int(float x)
{
	DoubleUnion u;
	u.D = x + DOUBLE2_FIX_MAGIC;

	return u.I[0];
}

inline unsigned long NextPow2(unsigned long X)
{
	unsigned long C = 1;
	while ((C < ULONG_MAX) && (C < X)) {
		C <<= 1;
	}
	return C;
}

// // this seems slower when profiling and only works for 32bit uint
//inline uint32_t NextPow2(uint32_t v)
//{
//	v--;
//	v |= v >> 1;
//	v |= v >> 2;
//	v |= v >> 4;
//	v |= v >> 8;
//	v |= v >> 16;
//	v++;
//	return v;
//
//}

inline bool IsPow2(unsigned long X)
{
	unsigned long C = 1;
	while ((C < ULONG_MAX) && (C < X)) {
		C <<= 1;
	}
	return (C == X);
}

// convert time*sampleRate to integer that is a multiple of 32
inline int ConvertTimeToSizeM32(const float time, const float sampleRate)
{
	int s = static_cast<int>(time*sampleRate) - 1; 
	return s - (s & 31) + 32;
}

inline int ConvertTimeToSizePow2(const float time, const float sampleRate)
{
	return static_cast<int>(NextPow2(static_cast<unsigned long>(time*sampleRate)));
}

inline int ConvertTimeToSizePow2Assert(const float time, const float sampleRate)
{
	assert(IsPow2(static_cast<unsigned long>(time*sampleRate)));
	return ConvertTimeToSizePow2(time, sampleRate);
}

inline float ConvertSizePow2ToTime(const int size, const float sampleRate)
{
	return static_cast<float>(NextPow2(static_cast<unsigned long>(size))) / sampleRate;
}

inline float ConvertSizePow2ToTimeAssert(const int size, const float sampleRate)
{
	assert(IsPow2(size));
	return ConvertSizePow2ToTime(size, sampleRate);
}

inline int Log2(int x)
{
	int y = 0;
	if (x <= 0) { return -1; }
	while (x >>= 1) { ++y; }
	return y;
}

inline float Pow2(int exponent)
{
	// Linux and Visual Studio 2015
#if defined(__linux__) || (_MSC_VER >= 1900)    
	DoubleUnion y = { 0 };
	// For speed we set the exponent directly in the double.
	y.I[1] = (1023 + exponent) << 20U;
	return static_cast<float>(y.D);
#else
	// On Windows pow(2.0, x) is really fast for Visual Studio 2012 
	return pow(2.0, exponent);
#endif
}

// 2^x with max relative error of 6.15%
inline float Pow2Coarse(float exponent)
{
	exponent = std::min(std::max(exponent, -31.f), 31.f);
	const int iexp = static_cast<int>(std::floor(exponent));
	float fexp = exponent + 1.0f - iexp;
	if (iexp<0) { fexp /= (1 << -iexp); }
	else { fexp *= (1 << iexp); }
	return fexp;
}

// 10^x with max relative error of 6.15%
inline float Pow10Coarse(float exponent)
{
	exponent *= 3.321928094887363f;
	return Pow2Coarse(exponent);
}

// e^x with max relative error of 6.15%
inline float PowECoarse(float exponent)
{
	exponent *= 1.442695040888963f;
	return Pow2Coarse(exponent);
}

// 10^(x/20) with max relative error of 6.15%
inline float Db2LinCoarse(float valueLin)
{
	valueLin *= 0.166096404744368f;
	return Pow2Coarse(valueLin);
}

inline int32_t Log2Int(double x)
{
	return static_cast<int32_t>((((reinterpret_cast<unsigned long long&>(x)) >> 52) & 0x7ff) - 1023);
}

// square root using second order taylor expansion on mantissa
inline float SqrtCoarse(float Value)
{
	int i;
	float x = std::frexp(Value, &i);
	x = x * (-0.198883056640625f * x + 0.881011962890625f) + 0.317169189453125f;
	return x * (1 << (i >> 1));
}

// Aproximation to tanh

inline float TanhApprox(const float x)
{
	const auto y = std::max(std::min(x, 3.f), -3.f);
	const auto y2 = y * y;
	return y * (27 + y2) / (27 + 9 * y2);
}

inline float TanhApprox2(const float x) 
{
	const auto y = std::max(std::min(x, 1.875f), -1.875f);
	const auto y2 = y * y;
	return y * (1.f - y2 * (0.189629629629630f - y2 * 0.016181728395062f)); // = y - 0.189629629629630f*y^3 + 0.016181728395062f*y^5
}

// Tan(x) approximation on interval [0,pi/4]. All other ranges are mirrors of this interval.
inline float TanApprox(const float x)
{
	return x / (1.003345165200669f - x * x * 0.352579063921458f);
}

inline float Clip1(const float x)
{
	return std::min(std::max(x, -1.f), 1.f);
}