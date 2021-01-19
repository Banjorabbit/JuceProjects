#pragma once
#pragma warning (push)
#pragma warning (disable : 4127 ) //4127 removes "unused variable" warnings from Eigen
#define EIGEN_DENSEBASE_PLUGIN "GetSizeOfEigen.h" //  member function added to Eigen DenseBase class to get dynamic memory size of array and matrices
#define EIGEN_MPL2_ONLY // don't allow LGPL licensed code from Eigen
#include <Eigen/Dense> // Eigen Library. Potentially redundant since this is also included in AsynchronousBase.h and InputOutput.h
#pragma warning (pop)
#include "Base.h" // Base class
#include "InputOutput.h" // Input/Output structs.
#include "Macros.h" // macros for member algorithms
#include "ConstrainedType.h" // ConstrainedType class

// AsynchronousBase inherits from Base and adds functions for asynchronous processing.
// It allows streaming data asynchronously (ie. with other buffer sizes than being used internally) 
// in and out of the algorithm, which is particularly useful for processing buffers from a 
// host that sends/receives buffers that:
// 1) are not supported by the algorithm (for instance not-power-of-2 buffer sizes).
// 2) produces sub-optimal results (for instance if the algorithm performs best at a certain buffer size).
// 3) are not known (or can change later) at the time of initialization.
//
// author: Kristian Timm Andersen
template<typename Talgo, typename Tinput = I::Real2D, typename Toutput = O::Real2D, typename Tpersistent = I::Real2D>
class AsynchronousBase : public Base<Talgo, Tinput, Toutput, Tpersistent>
{
	friend Base<Talgo, Tinput, Toutput, Tpersistent>;

public:

	// The variable AsynchronousBufferType AsynchronousBuffer must be defined in the coefficients of a class that derives from AsynchronousBase. 
	// This variable is checked during InitializeAsynchronous() to determine if the buffer size can be changed to match the external bufferSize. 
	// If it is equal to VARIABLE_SIZE, then it can change the buffer size during initialization. This will increase the chance that the algorithm can run synchronously.
	// If it is equal to FIXED_SIZE, then it will not change the buffer size.
	enum AsynchronousBufferType { FIXED_SIZE, VARIABLE_SIZE};

	// getters
	bool GetInitializedAsynchronous() const { return InitializedAsynchronous; }
	
	int GetLatencySamplesAsynchronous() const
	{ 
		if (InitializedAsynchronous) { return LatencyTotalSamples; }
		return static_cast<Talgo const&>(*this).GetLatencySamples();
	}

	// override this in derived class if it has internal latency
	int GetLatencySamples() const { return 0; } 
	
	int GetBufferSizeExternal() const
	{
		if (InitializedAsynchronous) { return BufferSizeExternal; }
		return static_cast<Talgo const&>(*this).GetCoefficients().BufferSize;
	}

	bool GetSynchronousProcessing() const
	{
		if (InitializedAsynchronous) { return SynchronousProcessing; }
		return true;
	}

	int GetNChannelsIn() const
	{
		return static_cast<Talgo const&>(*this).GetCoefficients().NChannelsIn;
	}

	// override this in derived class if NChannelsOut != NChannelsIn
	int GetNChannelsOut() const 
	{
		return static_cast<Talgo const&>(*this).GetCoefficients().NChannelsIn;
	}

	int GetBufferSize() const 
	{ 
		return static_cast<Talgo const&>(*this).GetCoefficients().BufferSize;
	}
	
	// override this in derived class if it shouldn't just accept the first suggested buffer size
	auto GetAsynchronousCoefficients(const std::vector<int> bufferSizesSuggested) const
	{
		auto c = static_cast<Talgo const&>(*this).GetCoefficients();
		c.BufferSize = bufferSizesSuggested[0];
		return c;
	}

	// overrides Base<Talgo>::Reset()
	void Reset() 
	{
		ResetAsynchronous();
		Base<Talgo, Tinput, Toutput, Tpersistent>::Reset();
	}

	// overrides Base<Talgo>::Initialize(c)
	template<typename Tcoefficients>
	auto Initialize(const Tcoefficients& c)
	{
		auto flag = Base<Talgo, Tinput, Toutput, Tpersistent>::Initialize(c);
		if (flag) 
		{ 
			InitializedAsynchronous = false; 
			BufferIn.resize(0, 0);
			BufferOut.resize(0, 0);
		}
		return flag;
	}

	// overrides Base<Talgo>::Initialize()
	auto Initialize()
	{
		auto flag = Base<Talgo, Tinput, Toutput, Tpersistent>::Initialize();
		if (flag) 
		{ 
			InitializedAsynchronous = false; 
			BufferIn.resize(0, 0);
			BufferOut.resize(0, 0);
		}
		return flag;
	}

	// overrides Base<Talgo>::GetAllocatedMemorySize()
	size_t GetAllocatedMemorySize() const
	{
		return Base<Talgo, Tinput, Toutput, Tpersistent>::GetAllocatedMemorySize() + BufferIn.GetAllocatedMemorySize() + BufferOut.GetAllocatedMemorySize();
	}

	template<typename Tcoefficients>
	bool InitializeAsynchronous(const Tcoefficients& c, const int bufferSizeExternal = 1)
	{
		auto& derived = static_cast<Talgo&>(*this);
		auto& derivedConst = static_cast<Talgo const&>(*this);

		BufferSizeExternal = bufferSizeExternal;
		auto cAsynchronous = c;
		if (c.AsynchronousBuffer == VARIABLE_SIZE)
		{
			derived.SetCoefficients(cAsynchronous);
			std::vector<int> bufferSizesSuggested = { bufferSizeExternal, std::max(bufferSizeExternal / 2, 1), std::max(bufferSizeExternal / 4, 1), std::max(bufferSizeExternal / 8, 1),
									 std::max(static_cast<int>(std::pow(2, std::round(std::log2(bufferSizeExternal)))), 1),
									 std::max(static_cast<int>(std::pow(2, std::ceil(std::log2(bufferSizeExternal)))), 1),
									 std::max(static_cast<int>(std::pow(2, std::floor(std::log2(bufferSizeExternal)))), 1) };
			cAsynchronous = derivedConst.GetAsynchronousCoefficients(bufferSizesSuggested);
		}
		
		bool flag = derived.Initialize(cAsynchronous);
		BufferSize = cAsynchronous.BufferSize;
		NChannelsIn = cAsynchronous.NChannelsIn;
		NChannelsOut = derivedConst.GetNChannelsOut();

		if (std::round(static_cast<float>(BufferSizeExternal) / BufferSize) == static_cast<float>(BufferSizeExternal) / BufferSize)
		{
			SynchronousProcessing = true;
		}
		else { SynchronousProcessing = false; }

		if (SynchronousProcessing) { LatencyTotalSamples = derivedConst.GetLatencySamples(); }
		else { LatencyTotalSamples = BufferSize + derivedConst.GetLatencySamples(); }
		
		BufferIn.resize(BufferSize, NChannelsIn);
		BufferOut.resize(BufferSize, NChannelsOut);

		ResetAsynchronous();

		if (flag) { InitializedAsynchronous = true; }

		return flag;
	}

	bool InitializeAsynchronous(const int bufferSizeExternal = 1)
	{ 
		auto c = static_cast<Talgo&>(*this).GetCoefficients();
		return InitializeAsynchronous(c, bufferSizeExternal);
	}
	
	// overrides Base<Talgo>::Process(xTime, yTime)
	void Process(Tinput xTime, Toutput yTime)
	{
		if (InitializedAsynchronous) { ProcessAsynchronous(xTime, yTime); }
		else { Base<Talgo, Tinput, Toutput, Tpersistent>::Process(xTime, yTime); }
	}

private:
	void ProcessAsynchronous(Tinput xTime, Toutput yTime)
	{
		if (SynchronousProcessing) { ProcessEntireBuffers(xTime, yTime); }
		else { ProcessSmallerBuffers(xTime, yTime); }
	}

	// Process the entire xTime buffer and put result in yTime. This function does not result in a delay
	void ProcessEntireBuffers(Tinput xTime, Toutput yTime)
	{
		int i = 0;
		// Process as many full buffers as possible
		for (i; i <= xTime.rows() - BufferSize; i += BufferSize)
		{
			Base<Talgo, Tinput, Toutput, Tpersistent>::Process(xTime.middleRows(i, BufferSize), yTime.middleRows(i, BufferSize));
		}
		// if we have been given a size that is not an integer multiple of BufferSizeInternal, zeropad and process. This is used if host does not give the expected length.
		int remainingSamples = std::min(static_cast<int>(xTime.rows()) - i, BufferSize);
		for (i; i < xTime.rows(); i += remainingSamples)
		{
			BufferIn.topRows(remainingSamples) = xTime.middleRows(i, remainingSamples);
			BufferIn.bottomRows(BufferSize - remainingSamples).setZero();
			Base<Talgo, Tinput, Toutput, Tpersistent>::Process(BufferIn, BufferOut);
			yTime.middleRows(i, remainingSamples) = BufferOut.topRows(remainingSamples);
			remainingSamples = std::min(static_cast<int>(xTime.rows()) - i, BufferSize);
		}
	}

	// Copy xTime buffer into BufferIn and process when it is full. While copying from xTime, also copy BufferOut (result of last Process) to yTime. This results in a delay equal to length of BufferSizeInternal.
	// This function is used when bufferSizeExpected is not an integer multiple of BufferSizeInternal (for instance if bufferSizeExpected is smaller than BufferSizeInternal, or 1.5x, 2.5x,... times larger).
	void ProcessSmallerBuffers(Tinput xTime, Toutput yTime)
	{
		for (auto i = 0; i < xTime.rows(); i++)
		{
			BufferIn.row(Index) = xTime.row(i);
			yTime.row(i) = BufferOut.row(Index);
			Index++;
			if (Index == BufferSize)
			{
				Base<Talgo, Tinput, Toutput, Tpersistent>::Process(BufferIn, BufferOut);
				Index = 0;
			}
		}
	}

	void ResetAsynchronous()
	{
		Index = 0;
		BufferIn.setZero();
		BufferOut.setZero();
	}

	typename I::GetType<Tinput>::type BufferIn; 
	typename O::GetType<Toutput>::type BufferOut;
	int LatencyTotalSamples = 0;
	int BufferSizeExternal = 0;
	int BufferSize = 0;
	int NChannelsIn = 0;
	int NChannelsOut = 0;
	bool SynchronousProcessing = true;
	bool InitializedAsynchronous = false;
	int Index = 0;
};