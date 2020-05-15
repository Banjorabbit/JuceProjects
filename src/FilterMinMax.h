#include "BaseClasses/PureCRTP.h"

struct O::FilterMinMax
{
	Real2D MinValue;
	Real2D MaxValue;
};


// StreamingMinMax finds the minimum and maximum value over a
// window with length C.Length for each new sample. It requires on
// average no more than 3 comparisons per sample. The algorithm uses 2
// double-ended queues for the minimum and maximum indices. A delay line
// is also used internally since in a true streaming application you need 
// to be able to call the algorithm succesively with new frames (or just 1 new sample),
// and you are not guaranteed that the input frame is as long as C.Length. 
// To be able to preallocate, the queues have been implemented as circular buffers. 
//
// A symmetric version of StreamingMinMax has been implemented below called "Filter". 
// It might be necessary to call the public ResetInitialValues function before Process, 
// if certain initial conditions are required. Also versions that only find Min/Max have been implemented.
//
// ref: Daniel Lemire, STREAMING MAXIMUM - MINIMUM FILTER USING NO MORE THAN THREE COMPARISONS PER ELEMENT
//
// author : Kristian Timm Andersen
class StreamingMinMax : public Base<StreamingMinMax, I::Real2D, O::FilterMinMax>
{
	friend Base<StreamingMinMax, I::Real2D, O::FilterMinMax>;

public:
	void ResetInitialValue(const float inputOld) { Reset(); D.InputOld.setConstant(inputOld); }
	void ResetInitialValue(const Eigen::Ref<const Eigen::ArrayXXf>& inputOld) { Reset(); D.InputOld = inputOld; }

private:
	struct Coefficients 
	{
		int Length = 10; // length of input to find min/max over
		int NChannels = 1;
	} C;

	struct Parameters {} P;
	struct Data 
	{
		Eigen::ArrayXXi MaxIndex, MinIndex;
		Eigen::ArrayXf InputOld;
		Eigen::ArrayXXf DelayLine;
		int Lf, Le, Uf, Ue;
		int Index;
		void Reset() // this Reset assumes previous values was 0. To select another value, use public ResetInitialValues(float inputOld)
		{
			InputOld.setZero();
			Lf = 0;
			Le = 0;
			Uf = 0;
			Ue = 0;
			Index = 0;
			MaxIndex.setZero();
			MinIndex.setZero();
			DelayLine.setZero();
		}
		bool InitializeMemory(const Coefficients& c)
		{
			InputOld.resize(c.NChannels);
			MaxIndex.resize(c.Length, c.NChannels);
			MinIndex.resize(c.Length, c.NChannels);
			DelayLine.resize(c.Length, c.NChannels);
			return true;
		}
		size_t GetAllocatedMemorySize() const 
		{
			size_t size = MinIndex.GetAllocatedMemorySize();
			size += MaxIndex.GetAllocatedMemorySize();
			size += InputOld.GetAllocatedMemorySize();
			size += DelayLine.GetAllocatedMemorySize();
			return size;
		}
		void OnParameterChange(const Parameters& p, const Coefficients& c) {}
	} D;

	void ProcessOn(Input x, Output y) 
	{
		for (auto channel = 0; channel < C.NChannels; channel++)
		{
			for (auto i = 0; i < x.rows(); i++)
			{
				D.DelayLine(D.Index, channel) = D.InputOld(channel);
				int indexNext = D.Index + 1;
				if (indexNext == C.Length) { indexNext = 0; }

				if (x(i, channel) > D.InputOld(channel))
				{
					D.MinIndex(D.Le, channel) = D.Index;
					
					D.Le += 1;
					if (D.Le == C.Length) { D.Le = 0; }
					if (indexNext == D.MinIndex(D.Lf, channel))
					{ 
						D.Lf += 1;
						if (D.Lf == C.Length) { D.Lf = 0; }
					}
					while (D.Uf != D.Ue)
					{
						auto Up = D.Ue - 1;
						if (Up < 0) { Up = C.Length - 1; }
						if (x(i, channel) <= D.DelayLine(D.MaxIndex(Up, channel), channel))
						{
							if (indexNext == D.MaxIndex(D.Uf, channel))
							{
								D.Uf += 1;
								if (D.Uf == C.Length) { D.Uf = 0; }
							}
							break;
						}
						D.Ue = Up;
					}
				}
				else 
				{
					D.MaxIndex(D.Ue, channel) = D.Index;
					D.Ue = (D.Ue + 1) % C.Length;
					if (indexNext == D.MaxIndex(D.Uf, channel))
					{
						D.Uf += 1;
						if (D.Uf == C.Length) { D.Uf = 0; }
					}
					while (D.Lf != D.Le)
					{
						auto Lp = D.Le - 1;
						if (Lp < 0) { Lp = C.Length - 1; }
						if (x(i, channel) >= D.DelayLine(D.MinIndex(Lp, channel), channel))
						{
							if (indexNext == D.MinIndex(D.Lf, channel))
							{
								D.Lf += 1;
								if (D.Lf == C.Length) { D.Lf = 0; }
							}
							break;
						}
						D.Le = Lp;
					}
				}
				if (D.Uf == D.Ue) { y.MaxValue(i, channel) = x(i, channel); }
				else { y.MaxValue(i, channel) = D.DelayLine(D.MaxIndex(D.Uf, channel), channel); }
				if (D.Lf == D.Le) { y.MinValue(i, channel) = x(i, channel); }
				else { y.MinValue(i, channel) = D.DelayLine(D.MinIndex(D.Lf, channel), channel); }
				D.InputOld(channel) = x(i, channel);
				D.Index = indexNext;
			}
		}
	}

	void ProcessOff(Input x, Output y) 
	{ 
		y.MaxValue.setZero(); 
		y.MinValue.setZero(); 
	}
};

// c.Length should be odd to be symmetric. If c.Length is even, minValue will correspond to c.Length+1, and maxValue will correspond to c.Length-1
class FilterMinMax : public Base<FilterMinMax, I::Real2D, O::FilterMinMax>
{
	friend Base<FilterMinMax, I::Real2D, O::FilterMinMax>;

public:
	StreamingMinMax Streaming;

private:
	struct Coefficients
	{
		int Length = 10; // length of input to find min/max over
		int NChannels = 2;
	} C;

	struct Parameters {} P;
	struct Data
	{
		int wHalf;
		void Reset() {	}
		bool InitializeMemory(const Coefficients& c)
		{
			wHalf = (c.Length - 1) / 2; 
			return true;
		}
		size_t GetAllocatedMemorySize() const
		{
			return 0;
		}
		void OnParameterChange(const Parameters& p, const Coefficients& c) {}
	} D;

	DEFINEMEMBERALGORITHMS(1, Streaming);

	auto InitializeMembers()
	{
		auto c = Streaming.GetCoefficients();
		c.Length = C.Length;
		c.NChannels = C.NChannels;
		return Streaming.Initialize(c);
	}

	void ProcessOn(Input x, Output y) 
	{
		Streaming.ResetInitialValue(x.row(0).transpose());
		// this output is discarded and is only used to update internal values
		Streaming.Process(x.topRows(D.wHalf), { y.MinValue.topRows(D.wHalf), y.MaxValue.topRows(D.wHalf) });
		// This is the shifted streaming filter, which creates a symmetric window
		Eigen::ArrayXXf xSymmetric(x.rows(), x.cols());
		xSymmetric.topRows(x.rows() - D.wHalf) = x.bottomRows(x.rows() - D.wHalf);
		xSymmetric.bottomRows(D.wHalf) = x.bottomRows<1>().replicate(D.wHalf, 1);
		Streaming.Process(xSymmetric, y);
	}

	void ProcessOff(Input x, Output y) 
	{ 
		y.MaxValue.setZero();
		y.MinValue.setZero();
	}
};

class StreamingMax : public Base<StreamingMax>
{
	friend Base<StreamingMax>;

public:
	void ResetInitialValue(const float inputOld) { Reset(); D.InputOld.setConstant(inputOld); }
	void ResetInitialValue(const Eigen::Ref<const Eigen::ArrayXXf>& inputOld) { Reset(); D.InputOld = inputOld; }

private:
	struct Coefficients
	{
		int Length = 10; // length of input to find min/max over
		int NChannels = 1;
	} C;

	struct Parameters {} P;
	struct Data
	{
		Eigen::ArrayXXf MaxValue;
		Eigen::ArrayXXi MaxIndex;
		Eigen::ArrayXf InputOld;
		int Uf, Ue;
		void Reset() // this Reset assumes previous values was 0. To select another value, use public ResetInitialValues(float inputOld)
		{
			InputOld.setZero();
			Uf = 0;
			Ue = 0;
			MaxValue.setZero();
			MaxIndex.setZero();
		}
		bool InitializeMemory(const Coefficients& c)
		{
			InputOld.resize(c.NChannels);
			MaxValue.resize(c.Length, c.NChannels);
			MaxIndex.resize(c.Length, c.NChannels);
			return true;
		}
		size_t GetAllocatedMemorySize() const
		{
			size_t size = MaxValue.GetAllocatedMemorySize();
			size += MaxIndex.GetAllocatedMemorySize();
			size += InputOld.GetAllocatedMemorySize();
			return size;
		}
		void OnParameterChange(const Parameters& p, const Coefficients& c) {}
	} D;

	void ProcessOn(Input x, Output y)
	{
		for (auto channel = 0; channel < C.NChannels; channel++)
		{
			for (auto i = 0; i < x.rows(); i++)
			{
				if (x(i, channel) > D.InputOld(channel))
				{
					while (D.Uf != D.Ue)
					{
						auto Up = D.Ue - 1;
						if (Up < 0) { Up = C.Length - 1; }
						if (x(i, channel) <= D.MaxValue(Up, channel))
						{
							if (i == C.Length + D.MaxIndex(D.Uf, channel))
							{
								D.Uf += 1;
								if (D.Uf == C.Length) { D.Uf = 0; }
							}
							break;
						}
						D.Ue = Up;
					}
				}
				else
				{
					D.MaxIndex(D.Ue, channel) = i - 1;
					D.MaxValue(D.Ue, channel) = D.InputOld(channel);
					D.Ue = (D.Ue + 1) % C.Length;
					if (i == C.Length + D.MaxIndex(D.Uf, channel))
					{
						D.Uf += 1;
						if (D.Uf == C.Length) { D.Uf = 0; }
					}
				}
				if (D.Uf == D.Ue) { y(i, channel) = x(i, channel); }
				else { y(i, channel) = D.MaxValue(D.Uf, channel); }
				D.InputOld(channel) = x(i, channel);
			}
		}
		// update Index to allow next Process() to continue as if they are successive frames
		D.MaxIndex -= static_cast<int>(x.rows());
	}

	void ProcessOff(Input x, Output y)
	{
		y.setZero();
	}
};

class StreamingMin : public Base<StreamingMin>
{
	friend Base<StreamingMin>;

public:
	void ResetInitialValue(const float inputOld) { Reset(); D.InputOld.setConstant(inputOld); }
	void ResetInitialValue(const Eigen::Ref<const Eigen::ArrayXXf>& inputOld) { Reset(); D.InputOld = inputOld; }

private:
	struct Coefficients
	{
		int Length = 10; // length of input to find min/max over
		int NChannels = 1;
	} C;

	struct Parameters {} P;
	struct Data
	{
		Eigen::ArrayXXf MinValue;
		Eigen::ArrayXXi MinIndex;
		Eigen::ArrayXf InputOld;
		int Lf, Le;
		void Reset() // this Reset assumes previous values was 0. To select another value, use public ResetInitialValues(float inputOld)
		{
			InputOld.setZero();
			Lf = 0;
			Le = 0;
			MinValue.setZero();
			MinIndex.setZero();
		}
		bool InitializeMemory(const Coefficients& c)
		{
			InputOld.resize(c.NChannels);
			MinValue.resize(c.Length, c.NChannels);
			MinIndex.resize(c.Length, c.NChannels);
			return true;
		}
		size_t GetAllocatedMemorySize() const
		{
			size_t size = MinValue.GetAllocatedMemorySize();
			size += MinIndex.GetAllocatedMemorySize();
			size += InputOld.GetAllocatedMemorySize();
			return size;
		}
		void OnParameterChange(const Parameters& p, const Coefficients& c) {}
	} D;

	void ProcessOn(Input x, Output y)
	{
		for (auto channel = 0; channel < C.NChannels; channel++)
		{
			for (auto i = 0; i < x.rows(); i++)
			{
				if (x(i, channel) > D.InputOld(channel))
				{
					D.MinIndex(D.Le, channel) = i - 1;
					D.MinValue(D.Le, channel) = D.InputOld(channel);
					D.Le += 1;
					if (D.Le == C.Length) { D.Le = 0; }
					if (i == C.Length + D.MinIndex(D.Lf, channel))
					{
						D.Lf += 1;
						if (D.Lf == C.Length) { D.Lf = 0; }
					}
				}
				else
				{
					while (D.Lf != D.Le)
					{
						auto Lp = D.Le - 1;
						if (Lp < 0) { Lp = C.Length - 1; }
						if (x(i, channel) >= D.MinValue(Lp, channel))
						{
							if (i == C.Length + D.MinIndex(D.Lf, channel))
							{
								D.Lf += 1;
								if (D.Lf == C.Length) { D.Lf = 0; }
							}
							break;
						}
						D.Le = Lp;
					}
				}
				if (D.Lf == D.Le) { y(i, channel) = x(i, channel); }
				else { y(i, channel) = D.MinValue(D.Lf, channel); }
				D.InputOld(channel) = x(i, channel);
			}
		}
		// update Index to allow next Process() to continue as if they are successive frames
		D.MinIndex -= static_cast<int>(x.rows());
	}

	void ProcessOff(Input x, Output y)
	{
		y.setZero();
	}
};

// c.Length should be odd to be symmetric. If c.Length is even, maxValue will correspond to c.Length-1
class FilterMax : public Base<FilterMax>
{
	friend Base<FilterMax>;

public:
	StreamingMax Streaming;

private:
	struct Coefficients
	{
		int Length = 10; // length of input to find min/max over
		int NChannels = 2;
	} C;

	struct Parameters {} P;
	struct Data
	{
		int wHalf;
		void Reset() {	}
		bool InitializeMemory(const Coefficients& c)
		{
			wHalf = (c.Length - 1) / 2;
			return true;
		}
		size_t GetAllocatedMemorySize() const
		{
			return 0;
		}
		void OnParameterChange(const Parameters& p, const Coefficients& c) {}
	} D;

	DEFINEMEMBERALGORITHMS(1, Streaming);

	auto InitializeMembers()
	{
		auto c = Streaming.GetCoefficients();
		c.Length = C.Length;
		c.NChannels = C.NChannels;
		return Streaming.Initialize(c);
	}

	void ProcessOn(Input x, Output y)
	{
		Streaming.ResetInitialValue(x.row(0).transpose());
		// this output is discarded and is only used to update internal values
		Streaming.Process(x.topRows(D.wHalf), y.topRows(D.wHalf) );
		// This is the shifted streaming filter, which creates a symmetric window
		Eigen::ArrayXXf xSymmetric(x.rows(), x.cols());
		xSymmetric.topRows(x.rows() - D.wHalf) = x.bottomRows(x.rows() - D.wHalf);
		xSymmetric.bottomRows(D.wHalf) = x.bottomRows<1>().replicate(D.wHalf, 1);
		Streaming.Process(xSymmetric, y);
	}

	void ProcessOff(Input x, Output y)
	{
		y.setZero();
	}
};

// c.Length should be odd to be symmetric. If c.Length is even, minValue will correspond to c.Length+1
class FilterMin : public Base<FilterMin>
{
	friend Base<FilterMin>;

public:
	StreamingMin Streaming;

private:
	struct Coefficients
	{
		int Length = 10; // length of input to find min/max over
		int NChannels = 2;
	} C;

	struct Parameters {} P;
	struct Data
	{
		int wHalf;
		void Reset() {	}
		bool InitializeMemory(const Coefficients& c)
		{
			wHalf = (c.Length - 1) / 2;
			return true;
		}
		size_t GetAllocatedMemorySize() const
		{
			return 0;
		}
		void OnParameterChange(const Parameters& p, const Coefficients& c) {}
	} D;

	DEFINEMEMBERALGORITHMS(1, Streaming);

	auto InitializeMembers()
	{
		auto c = Streaming.GetCoefficients();
		c.Length = C.Length;
		c.NChannels = C.NChannels;
		return Streaming.Initialize(c);
	}

	void ProcessOn(Input x, Output y)
	{
		Streaming.ResetInitialValue(x.row(0).transpose());
		// this output is discarded and is only used to update internal values
		Streaming.Process(x.topRows(D.wHalf), y.topRows(D.wHalf));
		// This is the shifted streaming filter, which creates a symmetric window
		Eigen::ArrayXXf xSymmetric(x.rows(), x.cols());
		xSymmetric.topRows(x.rows() - D.wHalf) = x.bottomRows(x.rows() - D.wHalf);
		xSymmetric.bottomRows(D.wHalf) = x.bottomRows<1>().replicate(D.wHalf, 1);
		Streaming.Process(xSymmetric, y);
	}

	void ProcessOff(Input x, Output y)
	{
		y.setZero();
	}
};