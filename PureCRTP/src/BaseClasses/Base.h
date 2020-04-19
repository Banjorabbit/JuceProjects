#pragma once
#include <vector> // used by VectorAlgo

// Base defines the interface for an algorithm using CRTP and exposes:
//	BASECLASSVERSION			- version number of base class
//  GetAllocatedMemorySize()	- return allocated memory when Initialize is called in bytes
//	GetStaticMemorySize()		- return static memory in bytes
//	Initialize(c)				- initializes the algorithm based on coefficients c
//	Initialize()				- initializes the algorithm with current constant values
//	GetCoefficients()				- return coefficients 
//	Reset()						- reset algorithm
//	GetParameters() 			- return parameters
//	SetParameters(p)			- set parameters in algorithm equal to p
//	Disable()					- disable algorithm
//	Enable()					- enable algorithm
//	Process(Input, Output)		- process input and write to output or bypass if disabled
//
// Template argument Talgo is the algorithm class and contains:
//	Coefficients					- struct containing all coefficients (changing these values requires calling Initialize(c))
//	Parameters					- struct containing all parameters (values that can be changed in real-time using SetParameters(p))
//	Data						- struct containing all data (internal values that can be reset and/or initialized)
//	Members						- struct containing all members (algorithms inside this algorithm)
//	ProcessOn(Input, Output)	- Process the input and write to output when enabled
//	ProcessOff(Input, Output)	- Bypass algorithm when disabled
//
// Template arguments Tinput and Toutput is the input/output types. See structs I and O for predefined types.
//
// The following is an empty template for a new algorithm called "AlgoCool" with input "I::Real2D" and output "O::Real2D":
//
// class AlgoCool : public Base<AlgoCool>
// {
//  friend Base<AlgoCool>;
//	struct Coefficients { } C;
//	struct Parameters { } P;
//	struct Data {
//		void Reset() { }
//		void InitializeMemory(const Coefficients& c) { }
//		size_t GetAllocatedMemorySize() const { return 0; }
//		void OnParameterChange(const Parameters& p, const Coefficients& c) { }
//	} D;
//	void ProcessOn(Input x, Output y) { }
//	void ProcessOff(Input x, Output y) { }
// };
//
// author: Kristian Timm Andersen

template<typename Talgo, typename Tinput = I::Real2D, typename Toutput = O::Real2D, typename Tpersistent = I::Real2D>
class Base
{
protected:
	struct Setup
	{
		Setup() = default;
		template<typename Tc, typename Tp>
		Setup(const Tc& c, const Tp& p) { Coefficients = c; Parameters = p; }
		typename Talgo::Coefficients Coefficients;
		typename Talgo::Parameters Parameters;
	};

public:
	
	static constexpr double PI = 3.14159265358979323846;
	static constexpr double BUTTERWORTH_Q = 0.707106781186548;
	static constexpr int BASECLASSVERSION = 1;

	size_t GetAllocatedMemorySize() const
	{
		auto& algo = static_cast<Talgo const&>(*this);
		return algo.D.GetAllocatedMemorySize() + algo.GetMemberAllocatedMemorySize();
	}

	static size_t GetStaticMemorySize() { return sizeof(Talgo); }

	template<typename Tcoefficients>
	bool Initialize(const Tcoefficients& c)
	{
		auto& algo = static_cast<Talgo&>(*this);
		algo.C = c;
		auto flag = algo.D.InitializeMemory(algo.C);
		algo.D.Reset();
		algo.D.OnParameterChange(algo.P, algo.C);
		flag &= algo.InitializeMembers();
		if (flag == true) { Initialized = true; return true; }
		else { Initialized = false; return false; }
	}

	bool Initialize() { return Initialize(static_cast<Talgo&>(*this).C); }
	
	bool Initialize(const Setup& s) 
	{
		SetParameters(s.Parameters); 
		return Initialize(s.Coefficients); 
	}

	void Reset()
	{
		auto& algo = static_cast<Talgo&>(*this);
		algo.D.Reset();
		algo.ResetMembers();
	}

	// get
	auto GetCoefficients() const { return static_cast<Talgo const&>(*this).C; } 
	auto GetParameters() const { return static_cast<Talgo const&>(*this).P; }
	auto GetSetup() const { return Setup(GetCoefficients(), GetParameters()); }
	auto GetCoefficientsAll() const { return static_cast<Talgo const&>(*this).GetCoefficientsAllImpl(); }
	auto GetParametersAll() const {	return static_cast<Talgo const&>(*this).GetParametersAllImpl();	}
	auto GetSetupAll() const { return static_cast<Talgo const&>(*this).GetSetupAllImpl(); }

	// set
	template<typename Tcoefficients>
	void SetCoefficients(const Tcoefficients& c)
	{
		static_cast<Talgo&>(*this).C = c;
		Initialized = false;
	}
	template<typename Tparameters>
	void SetParameters(const Tparameters& p)
	{
		auto& algo = static_cast<Talgo&>(*this);
		algo.P = p;
		algo.D.OnParameterChange(algo.P, algo.C);
	}
	template<typename Tsetup>
	void SetSetup(const Tsetup& setup) { SetCoefficients(setup.Coefficients); SetParameters(setup.Parameters); }
	template<typename Tcoefficients>
	void SetCoefficientsAll(const Tcoefficients& c) { static_cast<Talgo&>(*this).SetCoefficientsAllImpl(c); }
	template<typename Tparameters>
	void SetParametersAll(const Tparameters& p) { static_cast<Talgo&>(*this).SetParametersAllImpl(p); }
	template<typename Tsetup>
	void SetSetupAll(const Tsetup& s) { static_cast<Talgo&>(*this).SetSetupAllImpl(s); }

	void Disable() { Enabled = false; }
	void Enable() { Enabled = true; }

	using Input = const Tinput&;
	using Output = Toutput;

	inline void Process(Input x, Output y)
	{
		if (Enabled && Initialized) { static_cast<Talgo&>(*this).ProcessOn(x, y); }
		else { static_cast<Talgo&>(*this).ProcessOff(x, y); }
	}

	using InputPersistent = const Tpersistent&;
	void SetPersistentInput(InputPersistent input)
	{
		if (Initialized) { static_cast<Talgo&>(*this).ProcessPersistentInput(input); }
	}

	auto GetInitialized() const { return Initialized; }
	auto GetEnabled() const { return Enabled; }

	Base() = default;
	~Base() = default; // default destructor
	Base(const Base&) = delete; // delete copy constructor
	Base(Base&&) = default; // default move constructor necessary for resizing std::vector<Base>
	Base& operator=(const Base&) = delete; // delete copy assigment
	Base& operator=(Base&&) = delete; // delete move assignment

protected:

	// these functions will be hidden if defined in derived Talgo
	bool InitializeMembers() { return true; }
	void ProcessPersistentInput(InputPersistent input) { (void) input; }

	// these functions will be hidden if macro DEFINEMEMBERALGORITHMS(N,...) is declared in derived Talgo
	void ResetMembers() {}
	void SetupMembers() {}
	size_t GetMemberAllocatedMemorySize() const { return 0; }
	auto GetCoefficientsAllImpl() const { return GetCoefficients(); }
	auto GetParametersAllImpl() const { return GetParameters(); }
	auto GetSetupAllImpl() const { return GetSetup(); }
	template<typename Tcoefficients>
	void SetCoefficientsAllImpl(const Tcoefficients& c) { SetCoefficients(c); }
	template<typename Tparameters>
	void SetParametersAllImpl(const Tparameters& p) { SetParameters(p); }
	template<typename Tsetup>
	void SetSetupAllImpl(const Tsetup& s) { SetSetup(s); }

	// allows derived of Talgo to see Data
	auto GetPrivateData() { return &static_cast<Talgo&>(*this).D; }

private:
	
	bool Enabled = true;
	bool Initialized = false;
};