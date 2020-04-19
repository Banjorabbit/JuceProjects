#pragma once
#include "BaseClasses/PureCRTP.h"
#include "DesignFIRNonParametric.h"
#include <array>

// Process one sample of a FIR minimum phase filter.
//
// author: Kristian Timm Andersen

class FIRMinPhase : public Base<FIRMinPhase>
{
	friend Base<FIRMinPhase>;

public:
	DesignFIRNonParametric FilterCalculation;

private:
	struct Coefficients 
	{
		int FilterSize = 512;
		int NChannels = 2;
		float SampleRate = 16e3;
	} C;

	struct Parameters 
	{
		std::array<float, 20> Frequencies = { 100.f, 124.f, 155.f, 192.f, 239.f, 298.f, 370.f, 460.f, 573.f, 712.f, 886.f, 1102.f, 1370.f, 1704.f, 2120.f, 2637.f, 3279.f, 4079.f, 5073.f, 6310.f }; // 20 is max number
		std::array<float, 20> GaindB = { 0 }; // 20 is max size
		ConstrainedType<int> NGains = { 20, 1, 20 }; // number of frequencies/Gains
	} P;

	struct Data 
	{
		bool UpdateFilter;
		Eigen::ArrayXf Filter;
		Eigen::ArrayXXf FrameIn;
		int Index;
		void Reset() 
		{ 
			UpdateFilter = true; 
			FrameIn.setZero(); 
			Filter.setZero();
			Index = 0; }
		bool InitializeMemory(const Coefficients& c)	
		{
			Filter.resize(c.FilterSize);
			FrameIn.resize(c.FilterSize, c.NChannels);
			return true;
		}
		size_t GetAllocatedMemorySize() const { return Filter.GetAllocatedMemorySize() + FrameIn.GetAllocatedMemorySize(); }
		void OnParameterChange(const Parameters& p, const Coefficients& c) 
		{ 
			UpdateFilter = true;	
		}
	} D;

	DEFINEMEMBERALGORITHMS(1, FilterCalculation);

	auto InitializeMembers() 
	{
		auto s = FilterCalculation.GetSetup();
		s.Coefficients.FilterSize = C.FilterSize;
		s.Coefficients.SampleRate = C.SampleRate;
		return FilterCalculation.Initialize(s);
	}

	void ProcessOn(Input xTime, Output yTime) 
	{
		if (D.UpdateFilter) 
		{
			Eigen::Map<Eigen::ArrayXf> Frequencies(&P.Frequencies[0], P.NGains);
			Eigen::Map<Eigen::ArrayXf> GaindB(&P.GaindB[0], P.NGains);
			FilterCalculation.Process({ Frequencies, GaindB }, D.Filter); 
		}

		D.FrameIn.row(D.Index) = xTime.transpose();

		yTime = (D.FrameIn.topRows(D.Index)*D.Filter.tail(D.Index).replicate(1,C.NChannels)).colwise().sum().transpose();
		yTime += (D.FrameIn.bottomRows(C.FilterSize - D.Index)*D.Filter.tail(C.FilterSize - D.Index).replicate(1, C.NChannels)).colwise().sum().transpose();

		D.Index--;
		if (D.Index < 0) { D.Index = C.FilterSize - 1; }
		
	}
	void ProcessOff(Input xTime, Output yTime) { yTime = xTime; }
};