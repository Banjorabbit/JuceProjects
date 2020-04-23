#pragma once
#include <vector> 

// Wrapper for std::vector that exposes the same interface as Base except Process.
// https://stackoverflow.com/questions/44098116/a-c11-wrapper-class-on-top-of-std-vector
//
// author: Kristian Timm Andersen

template<typename Talgo>
class VectorAlgo
{
public:

	static constexpr int BASECLASSVERSION = Talgo::BASECLASSVERSION;

	size_t GetAllocatedMemorySize() const
	{
		size_t size = 0;
		for (auto& element : Vec)
		{
			size += element.GetAllocatedMemorySize();
		}
		return size;
	}

	static size_t GetStaticMemorySize() { return sizeof(VectorAlgo) + sizeof(Talgo)*Vec.size(); }

	template<typename Tcoefficients>
	bool Initialize(const Tcoefficients& c)
	{
		bool flag = true;
		for (auto& element : Vec) { flag &= element.Initialize(c); }
		return flag;
	}

	bool Initialize()
	{
		bool flag = true;
		for (auto& element : Vec) { flag &= element.Initialize(); }
		return flag;
	}

	void Reset() { for (auto& element : Vec) { element.Reset(); } }

	// get
	auto GetCoefficients() const
	{
		std::vector<decltype(Vec[0].GetCoefficients())> cVec(Vec.size());
		for (auto i = 0; i < Vec.size(); i++) { cVec[i] = Vec[i].GetCoefficients(); }
		return cVec;
	}
	auto GetParameters() const
	{
		std::vector<decltype(Vec[0].GetParameters())> pVec(Vec.size());
		for (auto i = 0; i < Vec.size(); i++) { pVec[i] = Vec[i].GetParameters(); }
		return pVec;
	}
	auto GetSetup() const
	{
		std::vector<decltype(Vec[0].GetSetup())> sVec(Vec.size());
		for (auto i = 0; i < Vec.size(); i++) { sVec[i] = Vec[i].GetSetup(); }
		return sVec;
	}

	auto GetCoefficientsAll() const
	{
		std::vector<decltype(Vec[0].GetCoefficientsAll())> cVec(Vec.size());
		for (auto i = 0; i < Vec.size(); i++) { cVec[i] = Vec[i].GetCoefficientsAll(); }
		return cVec;
	}

	auto GetParametersAll() const
	{
		std::vector<decltype(Vec[0].GetParametersAll())> pVec(Vec.size());
		for (auto i = 0; i < Vec.size(); i++) { pVec[i] = Vec[i].GetParametersAll(); }
		return pVec;
	}

	auto GetSetupAll() const
	{
		std::vector<decltype(Vec[0].GetSetupAll())> sVec(Vec.size());
		for (auto i = 0; i < Vec.size(); i++) { sVec[i] = Vec[i].GetSetupAll(); }
		return sVec;
	}

	// set
	template<typename Tcoefficients>
	void SetCoefficients(const Tcoefficients& c) { for (auto i = 0;i < size();i++) { Vec[i].SetCoefficients(c[i]); } }

	template<typename Tparameters>
	void SetParameters(const Tparameters& p) { for (auto i = 0;i < size();i++) { Vec[i].SetParameters(p[i]); } }

	template<typename Tsetup>
	void SetSetup(const Tsetup& setup) { for (auto i = 0;i < size();i++) { Vec[i].SetSetup(setup[i]); } }

	void Disable() { for (auto& element : Vec) { element.Disable(); } }
	void Enable() { for (auto& element : Vec) { element.Enable(); } }

	void resize(const size_t i) { Vec.resize(i); }
	size_t size() const { return Vec.size(); }
	Talgo& operator[](const size_t i) { return Vec[i]; }
	auto begin() { return Vec.begin(); }
	auto begin() const { return Vec.begin(); }
	auto end() { return Vec.end(); }
	auto end()   const { return Vec.end(); }

	auto GetInitialized() const
	{
		std::vector<decltype(Talgo::GetInitialized())> iVec(Vec.size());
		for (auto i = 0; i < Vec.size(); i++) { iVec[i] = Vec[i].Initialized(); }
		return iVec;
	}

	auto GetEnabled() const
	{
		std::vector<decltype(Talgo::GetEnabled())> gVec(Vec.size());
		for (auto i = 0; i < Vec.size(); i++) { gVec[i] = Vec[i].GetEnabled(); }
		return gVec;
	}

private:
	std::vector<Talgo> Vec;
};

