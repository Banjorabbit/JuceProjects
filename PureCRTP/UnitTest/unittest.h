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

namespace InterfaceTests // this namespace contains InterfaceTest and InterfaceAsynchronousTest and should be used in any unit test of an algorithm that derives from Base
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
		double duration = 0;
		int checkSum = 0;
		for (auto i = 0; i < 100;i++)
		{
			auto start = std::chrono::steady_clock::now();
			algo.Process(input, output);
			auto end = std::chrono::steady_clock::now();
			auto time = std::chrono::duration<double, std::micro>(end - start).count();
			duration += time / 100.f;//std::min(duration, time);
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
	
	template<typename Talgo>
	bool InitializeAsynchronousTest(typename Talgo::InputPersistent inputPersistent)
	{
		Talgo algo;
		auto c = algo.GetCoefficients();
		outputLog << "InitializeAsynchronousTest: \n";
		auto size = algo.GetStaticMemorySize();
		outputLog << "Static memory size: " << size << " bytes.\n";
		bool flag = algo.Initialize(c);
		if (!flag) { outputLog << "InitializeAsynchronousTest failed: Initialize(); returned false. \n\n"; return false; }
		auto sizeInit = algo.GetAllocatedMemorySize();
		outputLog << "Allocated memory size: " << sizeInit << " bytes.\n";
		flag = algo.InitializeAsynchronous(c, c.BufferSize);
		if (!flag) { outputLog << "InitializeAsynchronousTest failed: InitializeAsynchronous(); returned false. \n\n"; return false; }
		auto sizeAsync = algo.GetAllocatedMemorySize();
		outputLog << "Allocated asynchronous memory size: " << sizeAsync << " bytes.\n";
		if (sizeAsync <= sizeInit) { outputLog << "InitializeAsynchronousTest failed: Allocated memory did not increase. \n\n"; return false; }
		flag = algo.Initialize(c);
		if (!flag) { outputLog << "InitializeAsynchronousTest failed: Initialize(); returned false. \n\n"; return false; }
		algo.SetPersistentInput(inputPersistent);
		auto sizeOut = algo.GetAllocatedMemorySize();
		if (sizeInit == sizeOut) { outputLog << "InitializeAsynchronousTest successful.\n\n"; return true; }
		else { outputLog << "InitializeAsynchronousTest failed.\n\n"; return false; }
	}

	template<typename Talgo>
	bool SynchronousAsynchronousTest(typename Talgo::InputPersistent inputPersistent)
	{
		outputLog << "SynchronousAsynchronousTest: \n";
		Talgo algo, algoAsynchronous;

		
		auto c = algo.GetCoefficients();
		
		// initialize asynchronously and find a bufferSize that works for asynchronous algo
		auto flag = algoAsynchronous.InitializeAsynchronous(c, c.BufferSize);
		c = algoAsynchronous.GetCoefficients();
		flag &= algoAsynchronous.InitializeAsynchronous(c, c.BufferSize);

		// match BufferSize in algo
		flag &= algo.Initialize(c);

		if ((algo.GetSynchronousProcessing() == false) | (algoAsynchronous.GetSynchronousProcessing() == false))
		{
			return false;
		}
		// define input/output
		I::GetType<Talgo::Input>::type input(c.BufferSize, c.NChannelsIn);
		O::GetType<Talgo::Output>::type output(c.BufferSize, algo.GetNChannelsOut());
		O::GetType<Talgo::Output>::type outputAsynchronous(c.BufferSize, algoAsynchronous.GetNChannelsOut());
		algoAsynchronous.SetPersistentInput(inputPersistent);
		algo.SetPersistentInput(inputPersistent);

		for (int i = 0; i < 100; i++)
		{
			input.setRandom();
			// process algo asynchronous
			algoAsynchronous.Process(input, outputAsynchronous);
			// process algo
			algo.Process(input, output);
		}
		
		float error = (outputAsynchronous - output).abs2().sum();

		outputLog << "Synchronous delay: (" << algo.GetLatencySamples() << ", " << algoAsynchronous.GetLatencySamples() << ") samples.\n";
		outputLog << "Synchronous internal Buffersize: (" << algo.GetBufferSize() << ", " << algoAsynchronous.GetBufferSize() << ") samples.\n";
		if (error == 0.f && flag) { outputLog << "SynchronousAsynchronousTest successful.\n\n"; return true; }
		else { outputLog << "SynchronousAsynchronousTest failed.\n\n"; return false; }
	}

	template<typename Talgo>
	bool DifferentSizesAsynchronousTest(typename Talgo::InputPersistent inputPersistent)
	{
		outputLog << "DifferentSizesAsynchronousTest: \n";
		Talgo algo;
		auto c = algo.GetCoefficients();
		outputLog << "Default BufferSize: " << c.BufferSize << ".\n";
		// define input/output
		I::GetType<Talgo::Input>::type input2(c.BufferSize * 2, c.NChannelsIn);
		I::GetType<Talgo::Input>::type input125(static_cast<int>(c.BufferSize * 1.25f), c.NChannelsIn);
		O::GetType<Talgo::Output>::type output2(c.BufferSize * 2, algo.GetNChannelsOut());
		O::GetType<Talgo::Output>::type output125(static_cast<int>(c.BufferSize * 1.25f), algo.GetNChannelsOut());
		input2.setRandom();
		input125.setRandom();

		// setup to expect double internal buffer size
		bool flag = algo.InitializeAsynchronous(c, c.BufferSize * 2);
		algo.SetPersistentInput(inputPersistent);
		if (!flag) { outputLog << "DifferentSizesAsynchronousTest failed: InitializAsynchronouse(c, c.BufferSize*2); returned false. \n\n"; return false; }
		auto size = algo.GetAllocatedMemorySize();
		double duration = 1e10;
		int checkSum = 0;
		for (auto i = 0; i < 100;i++)
		{
			auto start = std::chrono::steady_clock::now();
			algo.Process(input2, output2);
			auto end = std::chrono::steady_clock::now();
			auto time = std::chrono::duration<double, std::micro>(end - start).count();
			duration = std::min(duration, time);
			void *readPtr = &output2;
			for (auto j = 0; j < sizeof(decltype(output2)); j++) { checkSum += static_cast<uint8_t*>(readPtr)[j]; } // checkSum ensures that all outputs are actually written
		}
		outputLog << "Execution time of processing 2x is: " << duration << "us.\n";
		outputLog << "Checksum: " << checkSum << ".\n";
		outputLog << "Delay: " << algo.GetLatencySamplesAsynchronous() << " samples.\n";
		outputLog << "Synchronous processing: " << algo.GetSynchronousProcessing() << ".\n";
		outputLog << "Internal Buffersize: " << algo.GetBufferSize() << " samples.\n";
		outputLog << "External Buffersize: " << algo.GetBufferSizeExternal() << " samples.\n";
		auto sizeOut = algo.GetAllocatedMemorySize();
		if (size != sizeOut) { outputLog << "DifferentSizesAsynchronousTest failed.\n\n"; return false; }

		// setup to expect 1.25x internal buffer size
		flag = algo.InitializeAsynchronous(c, static_cast<int>(c.BufferSize * 1.25f));
		algo.SetPersistentInput(inputPersistent);
		if (!flag) { outputLog << "DifferentSizesAsynchronousTest failed: InitializeAsynchronous(c, c.BufferSize * 1.25); returned false. \n\n"; return false; }
		size = algo.GetAllocatedMemorySize();
		duration = 1e10;
		checkSum = 0;
		for (auto i = 0; i < 100;i++)
		{
			auto start = std::chrono::steady_clock::now();
			algo.Process(input125, output125);
			auto end = std::chrono::steady_clock::now();
			auto time = std::chrono::duration<double, std::micro>(end - start).count();
			duration = std::min(duration, time);
			void *readPtr = &output125;
			for (auto j = 0; j < sizeof(decltype(output125)); j++) { checkSum += static_cast<uint8_t*>(readPtr)[j]; } // checkSum ensures that all outputs are actually written
		}
		outputLog << "Minimum execution time of processing 1.25x is: " << duration << "us.\n";
		outputLog << "Checksum: " << checkSum << ".\n";
		outputLog << "Delay: " << algo.GetLatencySamplesAsynchronous() << " samples.\n";
		outputLog << "Synchronous processing: " << algo.GetSynchronousProcessing() << ".\n";
		outputLog << "Internal Buffersize: " << algo.GetBufferSize() << " samples.\n";
		outputLog << "External Buffersize: " << algo.GetBufferSizeExternal() << " samples.\n";
		sizeOut = algo.GetAllocatedMemorySize();

		if (size == sizeOut) { outputLog << "DifferentSizesAsynchronousTest successful.\n\n"; return true; }
		else { outputLog << "DifferentSizesAsynchronousTest failed.\n\n"; return false; }
	}

	template<typename Talgo>
	bool AlgorithmInterfaceAsynchronousTest()
	{
		I::GetType<Talgo::InputPersistent>::type inputPersistent; // this variable is an empty Eigen array that is only used to match the number of arguments.
		return AlgorithmInterfaceAsynchronousTest<Talgo>(inputPersistent);
	}

	template<typename Talgo>
	bool AlgorithmInterfaceAsynchronousTest(typename Talgo::InputPersistent inputPersistent)
	{
		Talgo algo;
		auto c = algo.GetCoefficients();
		I::GetType<Talgo::Input>::type input(c.BufferSize, c.NChannelsIn);
		O::GetType<Talgo::Output>::type output(c.BufferSize, algo.GetNChannelsOut());
		input.setRandom();
		outputLog << "\n--------------------------------------------------------------------------------------\n";
		outputLog << "AlgorithmInterfaceTest: \n";
		auto successFlag = AlgorithmInterfaceTest<Talgo>(input, output, inputPersistent);
		outputLog << "\n";
		successFlag &= InitializeAsynchronousTest<Talgo>(inputPersistent);
		successFlag &= SynchronousAsynchronousTest<Talgo>(inputPersistent);
		successFlag &= DifferentSizesAsynchronousTest<Talgo>(inputPersistent);
		if (successFlag)
		{
			outputLog << "AlgorithmInterfaceTest successful." << std::endl;
		}
		else
		{
			outputLog << "AlgorithmInterfaceTest failed." << std::endl;
		}
		outputLog << "\n**************************************************************************************\n" << std::endl;
		return successFlag;
	}
}