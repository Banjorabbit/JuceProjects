#pragma once
#pragma warning (push)
#pragma warning (disable : 4127 ) //4127 removes "unused variable" warnings from Eigen
#define EIGEN_MPL2_ONLY // don't allow LGPL licensed code from Eigen
#include <Eigen/Dense> // Eigen Library.
#pragma warning (pop)

// author: Kristian Timm Andersen

// Eigen inputs should be const Eigen::Ref<const T>& types: 
struct I
{
	// define In type
	template<typename T>
	using In = const Eigen::Ref<const T>&;

	template<typename T>
	using InArray = In < Eigen::Array<T, Eigen::Dynamic, 1> >;

	template<typename T>
	using InArray2D = In < Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> >;

	using Real = In<Eigen::ArrayXf>;
	using Complex = In<Eigen::ArrayXcf>;
	using Bool = In<Eigen::Array<bool, Eigen::Dynamic, 1>>;
	using Real2D = In<Eigen::ArrayXXf>;
	using Complex2D = In<Eigen::ArrayXXcf>;
	using Bool2D = In<Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>>;
	using Boolean = const bool&;
	using Void = void*;

	using Real4X = In<Eigen::Array4Xf>;
	using Real4 = In<Eigen::Array4f>;
	using RealX2 = In<Eigen::ArrayX2f>;

	struct RealComplex { Real R; Complex C; };
	struct RealReal { Real R1; Real R2; };
	struct ComplexComplex { Complex C1; Complex C2; };

	struct BeamformerAdaptive;
	struct Deconvolver;
	struct EchoSuppressionCovariance;
	struct EchoCancellerMomentum;
	struct EchoCancellerNLMS;
	struct EchoCancellerNLMSMinPow;
	struct EchoCancellerToeplitz;
	struct CubicSpline;
	struct DesignFIRNonParametric;
	struct DesignIIRNonParametric;
	struct InterpolationCubic;
	struct InterpolationTemporal;
	struct InterpolateTonal;
	struct ToeplitzSolver;

	struct NonparametricEqualizerPersistent;

	// extract type using partial template specialization: https://stackoverflow.com/questions/301203/extract-c-template-parameters
	template<typename>
	struct GetType;

	template<typename T>
	struct GetType<In<T>> { typedef T type; };

	template<typename T>
	struct GetType<Eigen::Array<T, Eigen::Dynamic, 1>> { typedef T type; };

	template<typename T>
	struct GetType<Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>> { typedef T type; };

	template<typename>
	struct GetScalarType;

	template<typename T>
	struct GetScalarType<In<Eigen::Array<T, Eigen::Dynamic, 1>>> { typedef T type; };

	template<typename T>
	struct GetScalarType<In<Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>>> { typedef T type; };
};

// Eigen outputs should be Eigen::Ref<T> types: 
struct O
{
	// define Out type
	template<typename T>
	using Out = Eigen::Ref<T>;

	template<typename T>
	using OutArray = Out< Eigen::Array<T, Eigen::Dynamic, 1> >;

	template<typename T>
	using OutArray2D = Out< Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> >;

	using Real = Out<Eigen::ArrayXf>;
	using Complex = Out<Eigen::ArrayXcf>;
	using Bool = Out<Eigen::Array<bool, Eigen::Dynamic, 1>>;
	using Real2D = Out<Eigen::ArrayXXf>;
	using Complex2D = Out<Eigen::ArrayXXcf>;
	using Bool2D = Out<Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>>;
	
	struct RealComplex { Real R; Complex C; };
	struct RealReal { Real R1; Real R2; };
	struct ComplexComplex { Complex C1; Complex C2; };

	using Boolean = bool&;
	using Float = float&;
	using Void = void*;

	struct NoiseEstimationSPP;
	struct StateVariableFilter;
	struct DesignIIRMinPhase;
	struct DesignIIRNonParametric;
	struct FilterMinMax;

	// extract type using partial template specialization: https://stackoverflow.com/questions/301203/extract-c-template-parameters
	template<typename>
	struct GetType;

	template<typename T>
	struct GetType<Out<T>> { typedef T type; };

	template<typename T>
	struct GetType<Eigen::Array<T, Eigen::Dynamic, 1>> { typedef T type; };

	template<typename T>
	struct GetType<Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>> { typedef T type; };

	template<typename>
	struct GetScalarType;

	template<typename T>
	struct GetScalarType<Out<Eigen::Array<T, Eigen::Dynamic, 1>>> { typedef T type; };

	template<typename T>
	struct GetScalarType<Out<Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>>> { typedef T type; };
};