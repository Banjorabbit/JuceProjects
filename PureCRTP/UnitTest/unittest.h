#pragma once
#include "CppUnitTest.h"
#include "LoggerOutput.h"
#include <chrono>
#include "../src/BaseClasses/PureCRTP.h"

// using namespace in header files is bad practice, but this header file should only be included in unittest .cpp files.
using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace Eigen;

namespace { // empty namespace creates a unique LoggerOutput in each unittest .cpp file to prevent linker errors
	LoggerOutput outputLog;
}

namespace InterfaceTests // this namespace contains InterfaceTest and InterfaceStreamingTest and should be used in any unit test of an algorithm that derives from Base
{
	template<typename Talgo>
	bool VersionBaseTest()
	{
		Talgo algo;
		outputLog << "Base class version is: " << algo.BASECLASSVERSION << ".\n";
		return true;
	}

	template<typename Talgo>
	bool ProcessOffTest(typename Talgo::Input input, typename Talgo::Output output, typename Talgo::InputPersistent persistent)
	{
		Talgo algo;
		auto size = algo.GetAllocatedMemorySize();
		algo.SetPersistentInput(persistent);
		algo.Disable();
		double duration = 1e10;
		for (auto i = 0; i < 100;i++)
		{
			auto start = std::chrono::steady_clock::now();
			algo.Process(input, output);
			auto end = std::chrono::steady_clock::now();
			auto time = std::chrono::duration<double, std::micro>(end - start).count();
			duration = std::min(duration, time);
		}
		outputLog << "Execution time of ProcessOff is: " << duration << "us.\n";
		auto sizeOut = algo.GetAllocatedMemorySize();
		if (size == sizeOut) { outputLog << "ProcessOffTest successful.\n"; return true; }
		else { outputLog << "ProcessOffTest failed.\n"; return false; }
	}

	template<typename Talgo>
	bool ProcessOnTest(typename Talgo::Input input, typename Talgo::Output output, typename Talgo::InputPersistent persistent)
	{
		Talgo algo;
		bool flag = algo.Initialize();
		if (!flag) { outputLog << "ProcessOnTest failed: Initialize() returned false. \n"; return false; }
		auto size = algo.GetAllocatedMemorySize();
		algo.SetPersistentInput(persistent);
		double duration = 1e10;
		int checkSum = 0;
		for (auto i = 0; i < 100;i++)
		{
			auto start = std::chrono::steady_clock::now();
			algo.Process(input, output);
			auto end = std::chrono::steady_clock::now();
			auto time = std::chrono::duration<double, std::micro>(end - start).count();
			duration = std::min(duration, time);
			void *readPtr = &output;
			for (auto j = 0; j < sizeof(Talgo::Output); j++) { checkSum += static_cast<uint8_t*>(readPtr)[j]; } // checkSum ensures that all outputs are actually written
		}
		outputLog << "Execution time of ProcessOn is: " << duration << "us.\n";
		outputLog << "Checksum: " << checkSum << ".\n";
		auto sizeOut = algo.GetAllocatedMemorySize();
		if (size == sizeOut) { outputLog << "ProcessOnTest successful.\n"; return true; }
		else { outputLog << "ProcessOnTest failed.\n"; return false; }
	}

	template<typename Talgo>
	bool ResetTest()
	{
		Talgo algo;
		bool flag = algo.Initialize();
		if (!flag) { outputLog << "ResetTest failed: Initialize() returned false. \n"; return false; }
		auto size = algo.GetAllocatedMemorySize();
		algo.Reset();
		auto sizeOut = algo.GetAllocatedMemorySize();
		if (size == sizeOut) { outputLog << "ResetTest successful.\n"; return true; }
		else { outputLog << "ResetTest failed.\n"; return false; }
	}

	template<typename Talgo>
	bool GetSetTest()
	{
		Talgo algo;
		bool flag = algo.Initialize();
		if (!flag) { outputLog << "GetSetTest failed: Initialize() returned false. \n"; return false; }
		auto size = algo.GetAllocatedMemorySize();

		auto p = algo.GetParameters();
		algo.SetParameters(p);
		auto pAll = algo.GetParametersAll();
		algo.SetParametersAll(pAll);
		auto sizeOut = algo.GetAllocatedMemorySize();
		if (size != sizeOut) { outputLog << "GetSetTest failed. Parameter test failed.\n"; return false; }

		auto c = algo.GetCoefficients();
		algo.SetCoefficients(c);
		auto cAll = algo.GetCoefficientsAll();
		algo.SetCoefficientsAll(cAll);
		sizeOut = algo.GetAllocatedMemorySize();
		if (size != sizeOut) { outputLog << "GetSetTest failed. Coefficient test failed.\n"; return false; }

		auto s = algo.GetSetup();
		algo.SetSetup(s);
		auto sAll = algo.GetSetupAll();
		algo.SetSetupAll(sAll);
		sizeOut = algo.GetAllocatedMemorySize();
		if (size != sizeOut) { outputLog << "GetSetTest failed. Setup test failed.\n"; return false; }

		outputLog << "GetSetTest successful.\n";
		return true;
	}

	template<typename Talgo>
	bool InitializeTest()
	{
		Talgo algo;
		auto size = algo.GetStaticMemorySize();
		outputLog << "Static memory size: " << size << " bytes.\n";
		bool flag = algo.Initialize();
		if (!flag) { outputLog << "InitializeTest failed: Initialize() returned false. \n"; return false; }
		auto sizeInit = algo.GetAllocatedMemorySize();
		outputLog << "Allocated memory size: " << sizeInit << " bytes.\n";
		auto c = algo.GetCoefficients();
		flag = algo.Initialize(c);
		if (!flag) { outputLog << "InitializeTest failed: Initialize(c) returned false. \n"; return false; }
		auto sizeOut = algo.GetAllocatedMemorySize();
		if (sizeInit == sizeOut) { outputLog << "InitializeTest successful.\n"; return true; }
		else { outputLog << "InitializeTest failed.\n"; return false; }
	}

	template<typename Talgo>
	bool AlgorithmInterfaceTest(typename Talgo::Input input, typename Talgo::Output output)
	{
		I::GetType<typename Talgo::InputPersistent>::type temp; // this variable is an empty Eigen array that is only used to match the number of arguments.
		return AlgorithmInterfaceTest<Talgo>(input, output, temp);
	}

	template<typename Talgo>
	bool AlgorithmInterfaceTest(typename Talgo::Input input, typename Talgo::Output output, typename Talgo::InputPersistent persistent)
	{
		auto successFlag = VersionBaseTest<Talgo>();
		successFlag &= ProcessOffTest<Talgo>(input, output, persistent);
		successFlag &= ProcessOnTest<Talgo>(input, output, persistent);
		successFlag &= ResetTest<Talgo>();
		successFlag &= GetSetTest<Talgo>();
		successFlag &= InitializeTest<Talgo>();
		if (successFlag)
		{
			outputLog << "AlgorithmInterfaceTest successful." << std::endl;
		}
		else
		{
			outputLog << "AlgorithmInterfaceTest failed." << std::endl;
		}
		return successFlag;
	}

	template<typename TalgoStreaming>
	bool InitializeAsynchronousStreamingTest(I::Real2D xTime, float sampleRate, typename decltype(TalgoStreaming::Algo)::InputPersistent xPersistent)
	{
		TalgoStreaming algoStreaming;
		auto size = algoStreaming.GetStaticMemorySize();
		outputLog << "Static memory size: " << size << " bytes.\n";
		bool flag = algoStreaming.Initialize(static_cast<int>(xTime.rows()), static_cast<int>(xTime.cols()), sampleRate);
		if (!flag) { outputLog << "InitializeAsynchronousStreamingTest failed: Initialize(input.rows(), input.cols(), sampleRate); returned false. \n\n"; return false; }
		auto sizeInit = algoStreaming.GetAllocatedMemorySize();
		outputLog << "Allocated memory size: " << sizeInit << " bytes.\n";
		auto c = algoStreaming.Algo.GetCoefficients();
		flag = algoStreaming.Initialize(static_cast<int>(xTime.rows()), static_cast<int>(xTime.cols()), sampleRate, c);
		algoStreaming.Algo.SetPersistentInput(xPersistent);
		if (!flag) { outputLog << "InitializeAsynchronousStreamingTest failed: Initialize(input.rows(), input.cols(), sampleRate, c) returned false. \n\n"; return false; }
		auto sizeOut = algoStreaming.GetAllocatedMemorySize();
		if (sizeInit == sizeOut) { outputLog << "InitializeAsynchronousStreamingTest successful.\n\n"; return true; }
		else { outputLog << "InitializeAsynchronousStreamingTest failed.\n\n"; return false; }
	}

	template<typename TalgoStreaming>
	bool DifferentSizesAsynchronousStreamingTest(I::Real2D xTime, O::Real2D yTime, float sampleRate, typename decltype(TalgoStreaming::Algo)::InputPersistent xPersistent)
	{
		TalgoStreaming algoStreaming;

		// setup to expect half size of actual size
		bool flag = algoStreaming.Initialize(static_cast<int>(xTime.rows() / 2), static_cast<int>(xTime.cols()), sampleRate);
		algoStreaming.Algo.SetPersistentInput(xPersistent);
		if (!flag) { outputLog << "AsynchronousStreamingDifferentSizesTest failed: Initialize(input.rows(), input.cols()/2, sampleRate); returned false. \n\n"; return false; }
		auto size = algoStreaming.GetAllocatedMemorySize();
		double duration = 1e10;
		int checkSum = 0;
		for (auto i = 0; i < 100;i++)
		{
			auto start = std::chrono::steady_clock::now();
			algoStreaming.Process(xTime, yTime);
			auto end = std::chrono::steady_clock::now();
			auto time = std::chrono::duration<double, std::micro>(end - start).count();
			duration = std::min(duration, time);
			void *readPtr = &yTime;
			for (auto j = 0; j < sizeof(decltype(yTime)); j++) { checkSum += static_cast<uint8_t*>(readPtr)[j]; } // checkSum ensures that all outputs are actually written
		}
		outputLog << "Execution time of processing 2x is: " << duration << "us.\n";
		outputLog << "Checksum: " << checkSum << ".\n";
		outputLog << "Delay: " << algoStreaming.GetLatencyTotalSamples() << " samples.\n";
		outputLog << "Synchronous processing: " << algoStreaming.GetSynchronousProcessing() << ".\n";
		outputLog << "Internal Buffersize: " << algoStreaming.GetBufferSizeInternal() << " samples.\n";
		auto sizeOut = algoStreaming.GetAllocatedMemorySize();
		if (size != sizeOut) { outputLog << "AsynchronousStreamingDifferentSizesTest failed.\n\n"; return false; }

		// setup to expect 1.25x size of actual size
		flag = algoStreaming.Initialize(static_cast<int>(xTime.rows() * 1.25f), static_cast<int>(xTime.cols()), sampleRate);
		algoStreaming.Algo.SetPersistentInput(xPersistent);
		if (!flag) { outputLog << "AsynchronousStreamingDifferentSizesTest failed: Initialize(input.rows(), input.cols() * 1.25f); returned false. \n\n"; return false; }
		size = algoStreaming.GetAllocatedMemorySize();
		duration = 1e10;
		checkSum = 0;
		for (auto i = 0; i < 100;i++)
		{
			auto start = std::chrono::steady_clock::now();
			algoStreaming.Process(xTime, yTime);
			auto end = std::chrono::steady_clock::now();
			auto time = std::chrono::duration<double, std::micro>(end - start).count();
			duration = std::min(duration, time);
			void *readPtr = &yTime;
			for (auto j = 0; j < sizeof(decltype(yTime)); j++) { checkSum += static_cast<uint8_t*>(readPtr)[j]; } // checkSum ensures that all outputs are actually written
		}
		outputLog << "Minimum execution time of processing 0.8x is: " << duration << "us.\n";
		outputLog << "Checksum: " << checkSum << ".\n";
		outputLog << "Delay: " << algoStreaming.GetLatencyTotalSamples() << " samples.\n";
		outputLog << "Synchronous processing: " << algoStreaming.GetSynchronousProcessing() << ".\n";
		outputLog << "Internal Buffersize: " << algoStreaming.GetBufferSizeInternal() << " samples.\n";
		sizeOut = algoStreaming.GetAllocatedMemorySize();

		if (size == sizeOut) { outputLog << "AsynchronousStreamingDifferentSizesTest successful.\n\n"; return true; }
		else { outputLog << "AsynchronousStreamingDifferentSizesTest failed.\n\n"; return false; }
	}

	template<typename TalgoStreaming>
	bool SynchronousStreamingTest(I::Real2D xTime, O::Real2D yTime, float sampleRate, typename decltype(TalgoStreaming::Algo)::InputPersistent xPersistent)
	{
		// find a bufferSize that works for algoStreaming
		TalgoStreaming algoStreaming;
		ArrayXXf input = xTime;
		ArrayXXf output = yTime;
		bool flag = algoStreaming.Initialize(static_cast<int>(input.rows()), static_cast<int>(input.cols()), sampleRate);
		auto c = algoStreaming.Algo.GetCoefficients(); // get updated coefficients
		auto bufferSizeNew = algoStreaming.GetBufferSizeInternal();
		input.resize(bufferSizeNew, input.cols());
		input.setRandom();
		output.resize(bufferSizeNew, output.cols());

		// process algoStreaming
		flag &= algoStreaming.Initialize(static_cast<int>(input.rows()), static_cast<int>(input.cols()), sampleRate, c);
		algoStreaming.Algo.SetPersistentInput(xPersistent);
		algoStreaming.Process(input, output);
		ArrayXXf outputStreaming = output;

		// process Algo
		decltype(TalgoStreaming::Algo) algo;
		flag &= algo.Initialize(c);
		algo.SetPersistentInput(xPersistent);
		algo.Process(input, output);

		float error = (outputStreaming - output).abs2().sum();

		outputLog << "Synchronous delay: " << algoStreaming.GetLatencyTotalSamples() << " samples.\n";
		outputLog << "Synchronous processing: " << algoStreaming.GetSynchronousProcessing() << ".\n";
		outputLog << "Synchronous internal Buffersize: " << algoStreaming.GetBufferSizeInternal() << " samples.\n";
		if (error == 0.f && flag) { outputLog << "SynchronousStreamingTest successful.\n\n"; return true; }
		else { outputLog << "SynchronousStreamingTest failed.\n\n"; return false; }
	}

	template<typename TalgoStreaming>
	bool AlgorithmInterfaceStreamingTest(I::Real2D xTime, O::Real2D yTime, float sampleRate)
	{
		I::GetType<typename decltype(TalgoStreaming::Algo)::InputPersistent>::type temp; // this variable is an empty Eigen array that is only used to match the number of arguments.
		return AlgorithmInterfaceStreamingTest<TalgoStreaming>(xTime, yTime, sampleRate, temp);
	}

	template<typename TalgoStreaming>
	bool AlgorithmInterfaceStreamingTest(I::Real2D xTime, O::Real2D yTime, float sampleRate, typename decltype(TalgoStreaming::Algo)::InputPersistent xPersistent)
	{
		auto successFlag = AlgorithmInterfaceTest<decltype(TalgoStreaming::Algo)>(xTime, yTime, xPersistent);
		outputLog << "\n";
		successFlag &= InitializeAsynchronousStreamingTest<TalgoStreaming>(xTime, sampleRate, xPersistent);
		successFlag &= SynchronousStreamingTest<TalgoStreaming>(xTime, yTime, sampleRate, xPersistent);
		successFlag &= DifferentSizesAsynchronousStreamingTest<TalgoStreaming>(xTime, yTime, sampleRate, xPersistent);
		if (successFlag)
		{
			outputLog << "AlgorithmInterfaceStreamingTest successful." << std::endl;
		}
		else
		{
			outputLog << "AlgorithmInterfaceStreamingTest failed." << std::endl;
		}
		return successFlag;
	}
}