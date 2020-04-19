#pragma once
#include "BaseClasses/PureCRTP.h"

class CircBuffer : public Base<CircBuffer>
{
	friend Base<CircBuffer>;

public:
	// get data using operator()
	auto operator() (const int index, const int channel) const { return GetDelayValues(index, channel); }
	auto operator() (const float index, const int channel) const { return GetDelayValues(index, channel); }
	auto operator() (const int index) const { return GetDelayValues(index); }
	auto operator() (const float index) const { return GetDelayValues(index); }
	//push data
	void Push(I::Real2D xTime) { PushValues(xTime); }
	void PopLIFO(O::Real2D yTime) { PopValuesLIFO(yTime); }
	void PopFIFO(O::Real2D yTime) { PopValuesFIFO(yTime); }

private:

	struct Coefficients 
	{
		int NChannels = 2;
		int DelayLength = 1600;
	} C;

	struct Parameters {} P;

	struct Data 
	{
		Eigen::ArrayXXf DelayLine;
		int IndexStart, IndexFIFO;
		void Reset() 
		{
			IndexStart = static_cast<int>(DelayLine.rows())-1;
			IndexFIFO = 0;
			DelayLine.setZero(); 			
		}
		bool InitializeMemory(const Coefficients& c)	
		{ 
			DelayLine.resize(c.DelayLength, c.NChannels);
			return true;
		}
		size_t GetAllocatedMemorySize() const { return DelayLine.GetAllocatedMemorySize(); }
		void OnParameterChange(const Parameters& p, const Coefficients& c) {}
	} D;

	inline void IncrementIndex() { D.IndexStart++; if (D.IndexStart >= C.DelayLength) { D.IndexStart = 0; } }
	inline void IncrementIndex(const int increment) { D.IndexStart += increment; if (D.IndexStart >= C.DelayLength) { D.IndexStart -= C.DelayLength; } }
	inline void DecrementIndex() { D.IndexStart--; if (D.IndexStart < 0) { D.IndexStart = C.DelayLength - 1; } }
	inline void DecrementIndex(const int decrement) { D.IndexStart -= decrement; if (D.IndexStart < 0) { D.IndexStart += C.DelayLength; } }
	inline void IncrementFIFO() { D.IndexFIFO++; if (D.IndexFIFO >= C.DelayLength) { D.IndexFIFO = 0; } }
	inline void IncrementFIFO(const int increment) { D.IndexFIFO += increment; if (D.IndexFIFO >= C.DelayLength) { D.IndexFIFO -= C.DelayLength; } }


	// NOTE: Alternative to ProcessOn. Currently not used.
	void ProcessOn2(Input xTime, Output yTime) 
	{ 
		for (auto sample = 0; sample < xTime.rows(); sample++)
		{
			IncrementIndex();
			// for loop has been profiled to be faster than .col()
			for (auto channel = 0; channel < xTime.cols(); channel++)
			{
				yTime(sample, channel) = D.DelayLine(D.IndexStart, channel);
				D.DelayLine(D.IndexStart, channel) = xTime(sample, channel);
			}
		}
	}

	void ProcessOn(Input xTime, Output yTime)
	{
		IncrementIndex();
		const auto rowsFirst = std::min(C.DelayLength - D.IndexStart, static_cast<int>(xTime.rows()));
		const auto rowsSecond = std::max(0, static_cast<int>(xTime.rows()) - C.DelayLength + D.IndexStart);
		yTime.topRows(rowsFirst) = D.DelayLine.middleRows(D.IndexStart, rowsFirst);
		yTime.bottomRows(rowsSecond) = D.DelayLine.topRows(rowsSecond);
		D.DelayLine.middleRows(D.IndexStart, rowsFirst) = xTime.topRows(rowsFirst);
		D.DelayLine.topRows(rowsSecond) = xTime.bottomRows(rowsSecond);
		IncrementIndex(static_cast<int>(xTime.rows()) - 1);
	}

	void ProcessOff(Input xTime, Output yTime) { yTime = xTime; }

	// get values
	inline float GetDelayValues(const int index, const int channel) const
	{
		// find internal index
		auto newIndex = D.IndexStart - index;
		if (newIndex < 0) { newIndex += C.DelayLength; }
		return D.DelayLine(newIndex, channel);
	}
	inline float GetDelayValues(const float index, const int channel) const
	{
		// find internal index
		auto newIndex = D.IndexStart - index;
		if (newIndex < 0) { newIndex += C.DelayLength; }
		auto iindex = static_cast<int>(newIndex);
		const auto remainder = newIndex - iindex;
		auto value = D.DelayLine(iindex, channel);
		// find next index
		iindex++;
		if (iindex == C.DelayLength) { iindex = 0; }
		// linear interpolation
		value += remainder * (D.DelayLine(iindex, channel) - value);
		return value;
	}
	inline Eigen::ArrayXf GetDelayValues(const int index) const 
	{ 
		// find internal index
		auto newIndex = D.IndexStart - index;
		if (newIndex < 0) { newIndex += C.DelayLength; }
		return D.DelayLine.row(newIndex).transpose(); 
	}
	inline Eigen::ArrayXf GetDelayValues(const float index) const // this function cannot return a Reference, since it is calculating an interpolated value
	{ 
		// find internal index
		auto newIndex = D.IndexStart - index;
		if (newIndex < 0) { newIndex += C.DelayLength; }
		auto iindex = static_cast<int>(newIndex);
		const auto remainder = newIndex - iindex;
		Eigen::ArrayXf value = D.DelayLine.row(iindex).transpose();
		// find next index
		iindex++;
		if (iindex == C.DelayLength) { iindex = 0; }
		// linear interpolation
		value += remainder * (D.DelayLine.row(iindex).transpose() - value);
		return value;
	}

	// push values
	inline void PushValues(I::Real xTime) 
	{ // for loop has been profiled to be faster than .row()
		IncrementIndex();
		for (auto channel = 0; channel < xTime.cols(); channel++)
		{
			D.DelayLine(D.IndexStart, channel) = xTime(channel);
		}
	}
	inline void PushValues(I::Real2D xTime)
	{
		IncrementIndex();
		const auto rowsFirst = std::min(C.DelayLength - D.IndexStart, static_cast<int>(xTime.rows()));
		const auto rowsSecond = std::max(0, static_cast<int>(xTime.rows()) - C.DelayLength + D.IndexStart);
		D.DelayLine.middleRows(D.IndexStart, rowsFirst) = xTime.topRows(rowsFirst);
		D.DelayLine.topRows(rowsSecond) = xTime.bottomRows(rowsSecond);
		IncrementIndex(static_cast<int>(xTime.rows()) - 1);
	}

	// pop values
	inline void PopValuesLIFO(O::Real yTime)
	{
		for (auto channel = 0; channel < yTime.rows(); channel++)
		{
			yTime(channel) = D.DelayLine(D.IndexStart, channel);
		}
		DecrementIndex();
	}
	inline void PopValuesLIFO(O::Real2D yTime)
	{
		const auto rowsSecond = std::min(D.IndexStart + 1, static_cast<int>(yTime.rows()));
		const auto rowsFirst = std::max(0, static_cast<int>(yTime.rows()) - D.IndexStart - 1);
		yTime.topRows(rowsFirst) = D.DelayLine.bottomRows(rowsFirst);
		yTime.bottomRows(rowsSecond) = D.DelayLine.middleRows(D.IndexStart-rowsSecond+1,rowsSecond);
		DecrementIndex(static_cast<int>(yTime.rows()));
	}
	inline void PopValuesFIFO(O::Real yTime)
	{
		for (auto channel = 0; channel < yTime.rows(); channel++)
		{
			yTime(channel) = D.DelayLine(D.IndexFIFO, channel);
		}
		IncrementFIFO();
	}
	inline void PopValuesFIFO(O::Real2D yTime)
	{
		const auto rowsFirst = std::min(static_cast<int>(yTime.rows()), C.DelayLength - D.IndexFIFO);
		const auto rowsSecond = std::max(0, static_cast<int>(yTime.rows()) - C.DelayLength + D.IndexFIFO);
		yTime.topRows(rowsFirst) = D.DelayLine.middleRows(D.IndexFIFO, rowsFirst);
		yTime.bottomRows(rowsSecond) = D.DelayLine.topRows(rowsSecond);
		IncrementFIFO(static_cast<int>(yTime.rows()));
	}
};