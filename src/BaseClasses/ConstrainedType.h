#include <limits>

// Define ConstrainedType, which defines a min/max limit for the templated type.
//
// Author: Kristian Timm Andersen

template<typename T>
class ConstrainedType
{
public:
	ConstrainedType(const T& value, const T& vMin, const T& vMax) : VMin(vMin), VMax(vMax) { *this = value; }; // constructor, uses assignment operator to check value
	ConstrainedType& operator = (const T& value) // assignment operator
	{
		// // this commented-out code can be used to set a break point if necessary
		//if (value != (std::min)(VMax, value))
		//{
		//	int dummy = 1;
		//}
		//if (value != (std::max)(VMin, value))
		//{
		//	int dummy = 1;
		//}
		assert(value == (std::min)(VMax, value)); // catch in debug
		assert(value == (std::max)(VMin, value)); // catch in debug
		Value = (std::max)((std::min)(value, VMax), VMin);
		return *this;
	}
	operator T() const { return Value; } // cast to T type (used when assigning to other variable / expression)

	T Max() const { return VMax; }
	T Min() const { return VMin; }
private:
	T Value, VMin, VMax;
};

// This is used to get Eigen to accept that ConstrainedType<int> can be used as an integer, for instance: 
// ConstrainedType<int> N = {1,0,10};
// Eigen::ArrayXf X(N);
namespace std
{
	template<typename T>
	class numeric_limits<ConstrainedType<T>>
	{
	public:
		static constexpr bool is_integer = false;
		static constexpr bool is_signed = false;
	};

	template<>
	class numeric_limits<ConstrainedType<int>>
	{
	public:
		static constexpr bool is_integer = true;
		static constexpr bool is_signed = true;
	};
};