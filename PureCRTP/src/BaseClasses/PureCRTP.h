#pragma once
#pragma warning (push)
#pragma warning (disable : 4127 ) //4127 removes "unused variable" warnings from Eigen
#define EIGEN_DENSEBASE_PLUGIN "GetSizeOfEigen.h" //  member function added to Eigen DenseBase class to get dynamic memory size of array and matrices
#define EIGEN_MPL2_ONLY // don't allow LGPL licensed code from Eigen
#include <Eigen/Dense> // Eigen Library. Potentially redundant since this is also included in AsynchronousStreaming.h and InputOutput.h
#pragma warning (pop)
#include "Base.h" // Base class
#include "VectorAlgo.h" // VectorAlgo class
#include "AsynchronousStreaming.h" // AsynchronousStreaming class
#include "InputOutput.h" // Input/Output structs. Potentially redundant since this is also included in AsynchronousStreaming.h 
#include "Macros.h" // macros for member algorithms
#include "ConstrainedType.h" // ConstrainedType class

template<typename Talgo, typename Tinput = I::Complex2D, typename Toutput = O::Complex2D, typename Tpersistent = I::Complex2D>
using BaseFrequencyDomain = Base<Talgo, Tinput, Toutput, Tpersistent>;

// This is the main header file that includes the necessary files for developing algorithms using the PureCRTP library.
//
// author: Kristian Timm Andersen