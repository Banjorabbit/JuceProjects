#pragma once

// Assign values to array with an offset
//
// author: Kristian Timm Andersen
template<typename T, size_t m, size_t n>
void AssignArray(T(&array)[m], const T(&values)[n], const int offset = 0) { for (size_t i = 0; i < std::min(n, m - offset); i++) { array[i + offset] = values[i]; } }

// assign values from one pointer to another using Eigens fixed-size methods
// destination = target
template<typename Tx, typename Ty>
void inline assignEigen(Tx* destination, Ty* target, int length)
{
	int i = 0;
	for (i; i < length-3; i += 4, target += 4, destination += 4)
	{
		Eigen::Map<Eigen::Array<Tx, 4, 1>>(destination, 4) = Eigen::Map<Eigen::Array<Ty, 4, 1>>(target, 4);
	}
	for (i; i < length; i++, target++, destination++) { *destination = *target; }
}

// assign values from circular array to another using Eigens fixed-size methods
// destination(0:length1) = target(offset:length1)
// destination(length1:length1+length2) = target(0:length2)
template<typename T>
void inline assignCircularEigen(T* destination, T* target, int offset, int length1, int length2)
{
	assignEigen(destination, target + offset, length1);
	assignEigen(destination + length1, target, length2);
}

// multiply and assign values from one pointer to another using Eigens fixed-size methods
// destination = target * scalar
template<typename T, typename Tx, typename Ty>
void inline assignMultScalarEigen(T* destination, Tx* target, Ty value, int length)
{
	int i = 0;
	for (i; i < length - 3; i += 4, target += 4, destination += 4)
	{
		Eigen::Map<Eigen::Array<T, 4, 1>>(destination, 4) = Eigen::Map<Eigen::Array<Tx, 4, 1>>(target, 4) * value;
	}
	const auto len = length - i;
	Eigen::Map<Eigen::Array<T, Eigen::Dynamic, 1>>(destination, len) = Eigen::Map<Eigen::Array<Tx, Eigen::Dynamic, 1>>(target, len) * value;
}

// multiply and assign values from 2 pointers to another using Eigens fixed-size methods
// destination = target * array
template<typename T, typename Tx, typename Ty>
void inline assignMultEigen(T* destination, Tx* target, Ty* value, int length)
{
	int i = 0;
	for (i; i < length - 3; i += 4, target += 4, destination += 4, value += 4)
	{
		Eigen::Map<Eigen::Array<T, 4, 1>>(destination, 4) = Eigen::Map<Eigen::Array<Tx, 4, 1>>(target, 4) * Eigen::Map<Eigen::Array<Ty, 4, 1>>(value, 4);
	}
	const auto len = length - i;
	Eigen::Map<Eigen::Array<T, Eigen::Dynamic, 1>>(destination, len) = Eigen::Map<Eigen::Array<Tx, Eigen::Dynamic, 1>>(target, len) * Eigen::Map<Eigen::Array<Ty, Eigen::Dynamic, 1>>(value, len);
}

// multiply and assign values from 2 pointers to another using Eigens fixed-size methods
// destination = target * array
template<typename T, typename Tx, typename Ty>
void inline assignMultConjEigen(T* destination, Tx* target, Ty* value, int length)
{
	int i = 0;
	for (i; i < length - 3; i += 4, target += 4, destination += 4, value += 4)
	{
		Eigen::Map<Eigen::Array<T, 4, 1>>(destination, 4) = Eigen::Map<Eigen::Array<Tx, 4, 1>>(target, 4) * Eigen::Map<Eigen::Array<Ty, 4, 1>>(value, 4).conjugate();
	}
	const auto len = length - i;
	Eigen::Map<Eigen::Array<T, Eigen::Dynamic, 1>>(destination, len) = Eigen::Map<Eigen::Array<Tx, Eigen::Dynamic, 1>>(target, len) * Eigen::Map<Eigen::Array<Ty, Eigen::Dynamic, 1>>(value, len).conjugate();
}

// multiply and assign values from 2 pointers to another using Eigens fixed-size methods
// destination = target * array
template<typename T, typename Tx, typename Ty>
void inline assignDivEigen(T* destination, Tx* target, Ty* value, int length)
{
	int i = 0;
	for (i; i < length - 3; i += 4, target += 4, destination += 4, value += 4)
	{
		Eigen::Map<Eigen::Array<T, 4, 1>>(destination, 4) = Eigen::Map<Eigen::Array<Tx, 4, 1>>(target, 4) / Eigen::Map<Eigen::Array<Ty, 4, 1>>(value, 4);
	}
	const auto len = length - i;
	Eigen::Map<Eigen::Array<T, Eigen::Dynamic, 1>>(destination, len) = Eigen::Map<Eigen::Array<Tx, Eigen::Dynamic, 1>>(target, len) / Eigen::Map<Eigen::Array<Ty, Eigen::Dynamic, 1>>(value, len);
}

// multiply and assign values from 2 pointers to another using Eigens fixed-size methods
// destination = target * array
template<typename T, typename Tx>
void inline assignInvEigen(T* destination, Tx* target,  int length)
{
	int i = 0;
	for (i; i < length - 3; i += 4, target += 4, destination += 4)
	{
		Eigen::Map<Eigen::Array<T, 4, 1>>(destination, 4) = Eigen::Array<T,4,1>::Ones() / Eigen::Map<Eigen::Array<Tx, 4, 1>>(target, 4);
	}
	const auto len = length - i;
	Eigen::Map<Eigen::Array<T, Eigen::Dynamic, 1>>(destination, len) = Eigen::Array<T, Eigen::Dynamic, 1>::Ones(len) / Eigen::Map<Eigen::Array<Tx, Eigen::Dynamic, 1>>(target, len);
}

// multiply and assign values from 2 pointers to another using Eigens fixed-size methods
// destination = target * array
template<typename T, typename Tx, typename Ty>
void inline assignSubEigen(T* destination, Tx* target, Ty* value, int length)
{
	int i = 0;
	for (i; i < length - 3; i += 4, target += 4, destination += 4, value += 4)
	{
		Eigen::Map<Eigen::Array<T, 4, 1>>(destination, 4) = Eigen::Map<Eigen::Array<Tx, 4, 1>>(target, 4) - Eigen::Map<Eigen::Array<Ty, 4, 1>>(value, 4);
	}
	const auto len = length - i;
	Eigen::Map<Eigen::Array<T, Eigen::Dynamic, 1>>(destination, len) = Eigen::Map<Eigen::Array<Tx, Eigen::Dynamic, 1>>(target, len) - Eigen::Map<Eigen::Array<Ty, Eigen::Dynamic, 1>>(value, len);
}

// multiply and assign values from circular array to another using Eigens fixed-size methods
// destination(0:length1) = target(offset:length1) * value(0:length1)
// destination(length1:length1+length2) = target(0:length2) * value(length1:length1+length2)
template<typename T, typename Tx, typename Ty>
void inline assignMultCircularEigen(T* destination, Tx* target, Ty* value, int offset, int length1, int length2)
{
	assignMultEigen(destination, target + offset, value, length1);
	assignMultEigen(destination + length1, target, value + length1, length2);
}


// Abs2 and assign values from a pointer to another using Eigens fixed-size methods
// destination = target * array
template<typename T, typename Tx>
void inline assignAbs2Eigen(T* destination, Tx* target, int length)
{
	int i = 0;
	for (i; i < length - 3; i += 4, target += 4, destination += 4)
	{
		Eigen::Map<Eigen::Array<T, 4, 1>>(destination, 4) = Eigen::Map<Eigen::Array<Tx, 4, 1>>(target, 4).abs2();
	}
	const auto len = length - i;
	Eigen::Map<Eigen::Array<T, Eigen::Dynamic, 1>>(destination, len) = Eigen::Map<Eigen::Array<Tx, Eigen::Dynamic, 1>>(target, len).abs2();
}

// Add and assign values from a pointer to another using Eigens fixed-size methods
// destination += targe
template<typename T, typename Tx>
void inline addAssignEigen(T* destination, Tx* target, int length)
{
	int i = 0;
	for (i; i < length - 3; i += 4, target += 4, destination += 4)
	{
		Eigen::Map<Eigen::Array<T, 4, 1>>(destination, 4) += Eigen::Map<Eigen::Array<Tx, 4, 1>>(target, 4);
	}
	const auto len = length - i;
	Eigen::Map<Eigen::Array<T, Eigen::Dynamic, 1>>(destination, len) += Eigen::Map<Eigen::Array<Tx, Eigen::Dynamic, 1>>(target, len);
}

// Add and assign values from a pointer to another using Eigens fixed-size methods
// destination += targe
template<typename T, typename Tx>
void inline addScalarAssignEigen(T* destination, Tx value, int length)
{
	int i = 0;
	for (i; i < length - 3; i += 4, destination += 4)
	{
		Eigen::Map<Eigen::Array<T, 4, 1>>(destination, 4) += value;
	}
	const auto len = length - i;
	Eigen::Map<Eigen::Array<T, Eigen::Dynamic, 1>>(destination, len) += value;
}

// Add and assign values from a pointer to another using Eigens fixed-size methods
// destination += target * array
template<typename T, typename Tx, typename Ty>
void inline addAssignMultEigen(T* destination, Tx* target, Ty* value, int length)
{
	int i = 0;
	for (i; i < length - 3; i += 4, target += 4, destination += 4, value += 4)
	{
		Eigen::Map<Eigen::Array<T, 4, 1>>(destination, 4) += Eigen::Map<Eigen::Array<Tx, 4, 1>>(target, 4) * Eigen::Map<Eigen::Array<Tx, 4, 1>>(value, 4);
	}
	const auto len = length - i;
	Eigen::Map<Eigen::Array<T, Eigen::Dynamic, 1>>(destination, len) += Eigen::Map<Eigen::Array<Tx, Eigen::Dynamic, 1>>(target, len) * Eigen::Map<Eigen::Array<Tx, Eigen::Dynamic, 1>>(value, len);
}

// Add and assign values from a pointer to another using Eigens fixed-size methods
// destination += target * array
template<typename T, typename Tx, typename Ty>
void inline addAssignMultConjEigen(T* destination, Tx* target, Ty* value, int length)
{
	int i = 0;
	for (i; i < length - 3; i += 4, target += 4, destination += 4, value += 4)
	{
		Eigen::Map<Eigen::Array<T, 4, 1>>(destination, 4) += Eigen::Map<Eigen::Array<Tx, 4, 1>>(target, 4) * Eigen::Map<Eigen::Array<Tx, 4, 1>>(value, 4).conjugate();
	}
	const auto len = length - i;
	Eigen::Map<Eigen::Array<T, Eigen::Dynamic, 1>>(destination, len) += Eigen::Map<Eigen::Array<Tx, Eigen::Dynamic, 1>>(target, len) * Eigen::Map<Eigen::Array<Tx, Eigen::Dynamic, 1>>(value, len).conjugate();
}

// Subtract and assign values from a pointer to another using Eigens fixed-size methods
// destination = target * array
template<typename T, typename Tx>
void inline subAssignEigen(T* destination, Tx* target, int length)
{
	int i = 0;
	for (i; i < length - 3; i += 4, target += 4, destination += 4)
	{
		Eigen::Map<Eigen::Array<T, 4, 1>>(destination, 4) -= Eigen::Map<Eigen::Array<Tx, 4, 1>>(target, 4);
	}
	const auto len = length - i;
	Eigen::Map<Eigen::Array<T, Eigen::Dynamic, 1>>(destination, len) -= Eigen::Map<Eigen::Array<Tx, Eigen::Dynamic, 1>>(target, len);
}

// Multiply and assign values from a pointer to another using Eigens fixed-size methods
// destination = target * array
template<typename T, typename Tx>
void inline multAssignEigen(T* destination, Tx* target, int length)
{
	int i = 0;
	for (i; i < length - 3; i += 4, target += 4, destination += 4)
	{
		Eigen::Map<Eigen::Array<T, 4, 1>>(destination, 4) *= Eigen::Map<Eigen::Array<Tx, 4, 1>>(target, 4);
	}
	const auto len = length - i;
	Eigen::Map<Eigen::Array<T, Eigen::Dynamic, 1>>(destination, len) *= Eigen::Map<Eigen::Array<Tx, Eigen::Dynamic, 1>>(target, len);
}

// Multiply scalar and assign values from a pointer to another using Eigens fixed-size methods
// destination = target * array
template<typename T, typename Tx>
void inline multScalarAssignEigen(T* destination, Tx value, int length)
{
	int i = 0;
	for (i; i < length - 3; i += 4, destination += 4)
	{
		Eigen::Map<Eigen::Array<T, 4, 1>>(destination, 4) *= value;
	}
	const auto len = length - i;
	Eigen::Map<Eigen::Array<T, Eigen::Dynamic, 1>>(destination, len) *= value;
}
