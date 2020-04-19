#include "CppUnitTest.h"

#include "../src/CircBuffer.h"
#include "../src/CubicSpline.h"
#include "../src/DesignIIRMinPhase.h"
#include "../src/DesignIIRNonParametric.h"
#include "../src/DesignFIRNonParametric.h"
#include "../src/FFT.h"
#include "../src/Filterbank.h"
#include "../src/FIRMinPhase.h"
#include "../src/InterpolationCubic.h"
#include "../src/Upsampling2XCubic.h"

#include "../src/FrequencyDomain/BeamformerAdaptive.h"
#include "../src/FrequencyDomain/EchoCancellerMomentum.h"
#include "../src/FrequencyDomain/EchoSuppressionCovariance.h"
#include "../src/FrequencyDomain/GainCalculation.h"
#include "../src/FrequencyDomain/InterpolationTemporal.h"
#include "../src/FrequencyDomain/MinPhaseSpectrum.h"
#include "../src/FrequencyDomain/NoiseEstimation.h"
#include "../src/FrequencyDomain/NoiseSuppression.h"
#include "../src/FrequencyDomain/TonalDetection.h"
#include "../src/FrequencyDomain/VoiceActivationDetection.h"

#include "../src/AsynchronousStreaming/PipelinePreProcessing.h"
#include "../src/AsynchronousStreaming/FFTConvolution.h"
#include "../src/AsynchronousStreaming/IIR2ndDF.h"
#include "../src/AsynchronousStreaming/IIR2ndNonLin.h"
#include "../src/AsynchronousStreaming/IIR2ndTimeVarying.h"
#include "../src/AsynchronousStreaming/NonparametricEqualizer.h"
#include "../src/AsynchronousStreaming/PitchShift.h"
#include "../src/AsynchronousStreaming/VirtualizationHeadphones.h"

#include "LoggerOutput.h"
#include "../Utilities/AudioFile.h"
#include "../Utilities/FormatHandling.h"
#include "../Utilities/pffft.h"

#include <chrono>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace Eigen;

LoggerOutput outputLog;
//
//namespace Tuple
//{
//
//	TEST_CLASS(AlgorithmTupleTest)
//	{
//	public:
//
//		TEST_METHOD(Initialization)
//		{
//			AlgorithmTuple<NoiseEstimationSPP, std::vector<NoiseEstimationSPP>, NoiseEstimationSPP> test;
//			test.Get<1>().resize(2);
//			test.Initialize();
//			auto c = test.GetCoefficients();
//			std::get<0>(c).FilterbankRate = 110;
//			std::get<1>(c)[0].FilterbankRate = 120;
//			std::get<1>(c)[1].FilterbankRate = 140;
//			std::get<2>(c).FilterbankRate = 160;
//			test.Initialize(c);
//
//			auto cNew = test.GetCoefficients();
//			Assert::IsTrue(std::get<0>(cNew).FilterbankRate == 110);
//			Assert::IsTrue(std::get<1>(cNew)[0].FilterbankRate == 120);
//			Assert::IsTrue(std::get<1>(cNew)[1].FilterbankRate == 140);
//			Assert::IsTrue(std::get<2>(cNew).FilterbankRate == 160);
//		}
//		TEST_METHOD(GetCoefAndParam)
//		{
//			AlgorithmTuple<NoiseEstimationSPP, std::vector<NoiseEstimationSPP>, NoiseEstimationSPP> test;
//			test.Get<1>().resize(2);
//			test.Initialize();
//			auto c = test.GetCoefficients();
//			auto c0 = std::get<0>(c);
//			auto c1 = std::get<1>(c);
//			auto c2 = std::get<2>(c);
//
//			Assert::IsTrue(c0.FilterbankRate == c1[1].FilterbankRate);
//			Assert::IsTrue(c2.FilterbankRate == c1[0].FilterbankRate);
//
//			auto p = test.GetParameters();
//			auto p0 = std::get<0>(p);
//			auto p1 = std::get<1>(p);
//			auto p2 = std::get<2>(p);
//
//			Assert::IsTrue(p2.ActivityMeanTConstant == p1[1].ActivityMeanTConstant);
//			Assert::IsTrue(p0.ActivityMeanTConstant == p1[0].ActivityMeanTConstant);
//		}
//
//		TEST_METHOD(GetAllocatedMemorySize)
//		{
//			AlgorithmTuple<NoiseEstimationSPP, std::vector<NoiseEstimationSPP>> test;
//			test.Reset();
//			auto size = test.GetAllocatedMemorySize();
//			test.Get<1>().resize(2);
//			auto size0 = test.GetAllocatedMemorySize();
//			test.Initialize();
//			auto size1 = test.GetAllocatedMemorySize();
//
//			Assert::IsTrue(size == 0);
//			Assert::IsTrue(size0 == sizeof(NoiseEstimationSPP) * 2);
//			Assert::IsTrue(size1 > size0);
//		}
//
//		TEST_METHOD(SizeOf)
//		{
//			AlgorithmTuple<> test1;
//			AlgorithmTuple<NoiseEstimationSPP> test2;
//			NoiseEstimationSPP test3;
//			Assert::IsTrue(sizeof(test1) == 1);
//			Assert::IsTrue(sizeof(test2) == sizeof(test3));
//		}
//
//		TEST_METHOD(SetParametersTest)
//		{
//			AlgorithmTuple<NoiseEstimationSPP> test;
//			test.Initialize();
//			auto p = test.GetParameters();
//			std::get<0>(p).ActivityMeanTConstant = 1.f;
//			test.SetParameters(p);
//			auto pNew = test.GetParameters();
//			Assert::IsTrue(std::get<0>(pNew).ActivityMeanTConstant == 1.f);
//		}
//
//		TEST_METHOD(GetCoefficientsAll)
//		{
//			VoiceActivationDetection ref;
//			auto cRef = ref.GetCoefficients();
//			auto cTest = ref.GetCoefficientsAll();
//			
//			//Assert::IsTrue(cVecNoiseEst.size() == cRef.NChannels);
//			//for (auto i = 0;i < cRef.NChannels;i++)
//			//{
//			//	Assert::IsTrue(cVecNoiseEst[i].FilterbankRate == cRef.FilterbankRate);
//			//	Assert::IsTrue(cVecNoiseEst[i].NBands == cRef.NBands);
//			//	Assert::IsTrue(cVecNoiseEst[i].NChannels == 1);
//			//}
//
//			//auto pTuple = test.GetParametersAll();
//		}
//
//		TEST_METHOD(SizeOfEmptyTuple)
//		{
//			AlgorithmTuple<> test;
//			AlgorithmTuple<NoiseEstimationSPP> test2;
//			auto flag = sizeof(test) == sizeof(test2);
//			outputLog << "Equal: " << flag << std::endl;
//		}
//	};
//}

namespace Eigen
{
	TEST_CLASS(EigenMethodsTest)
	{
		TEST_METHOD(SizesTest)
		{
			ArrayXXf A;
			outputLog << "Size of A: " << A.size() << std::endl;
		}

		TEST_METHOD(ArrayXfTest)
		{
			ArrayXf test(4);
			outputLog << "Columns: " << test.cols() << std::endl;
		}
		TEST_METHOD(RefTest)
		{
			ArrayXXf test;
			ArrayXXcf testComplex;
			Ref<ArrayXXf> testRef = test;
			//const Ref<const ArrayXXf>& testRef2 = Map< const Ref<const ArrayXXf>,0>(&testComplex.real()(0,0),0,0);
			//decltype(testRef2.const_cast_derived()) test2 = test;
			outputLog << testRef.const_cast_derived() << std::endl;
		}
	};

	TEST_CLASS(FFTComparisons)
	{
		TEST_METHOD(PFFFT)
		{
			pffft_transform_t datatypeFFT = PFFFT_REAL;
			pffft_direction_t directionFFT = PFFFT_FORWARD;
			int N = 32;
			PFFFT_Setup *setup = pffft_new_setup(N, datatypeFFT);
			ArrayXf input(N);
			input = ArrayXf::LinSpaced(N, 0.f, N - 1.f);
			/*input.head(8) << 1, 2, 3, 2, 1, 2, 3, 2;
			input.segment(8,8) << 1, 2, 3, 2, 1, 2, 3, 2;
			input.segment(16, 8) << 1, 2, 3, 2, 1, 2, 3, 2;
			input.segment(24, 8) << 1, 2, 3, 2, 1, 2, 3, 2;*/
			outputLog << input << std::endl;
			ArrayXcf output(N/2+1);
			
			pffft_transform_ordered(setup,&input(0),&input(0), NULL, directionFFT);
			output(0) = input(0);
			output.tail(1) = input(1);
			std::memcpy(&output.real()(1), &input(2), (N - 2) * sizeof(float));
			outputLog << output << std::endl;

		}
		
	};

	TEST_CLASS(Efficiency)
	{
#if NDEBUG

		TEST_METHOD(EigenValues)
		{
			srand((unsigned int)time(0));
			EigenSolver<MatrixXf> es;
			int N = 20;
			es = EigenSolver<MatrixXf>::EigenSolver(N);

			MatrixXf mat(N,N);
			mat.topRows<1>().setRandom();
			mat.block(1, 0, N - 1, N - 1).setIdentity();
			mat.block(1, N - 1, N - 1, 1).setZero();

			es.compute(mat);

			double duration = 1e10;
			int checkSum = 0;
			for (auto i = 0; i < 100;i++)
			{
				auto start = std::chrono::steady_clock::now();
				es.compute(mat);
				auto end = std::chrono::steady_clock::now();
				auto time = std::chrono::duration<double, std::micro>(end - start).count();
				duration = std::min(duration, time);
			}
			outputLog << "Execution time of ProcessOn is: " << duration << "us.\n";
			MatrixXcf ev = es.eigenvalues();
			outputLog << "Matrix: " << mat.topRows<1>() << std::endl;
			outputLog << "Eigen values: " << ev.transpose() << std::endl;

			for (auto i = 0; i < ev.size();i++)
			{
				if (ev(i).imag() == 0) { outputLog << i << ", "; }
			}
			outputLog << std::endl;
			//ArrayX2f B(static_cast<int>((N + 1) / 2));
			//int index = 0;
			//for (auto i = 0;i < B.size();i++)
			//{

			//}
		}
		TEST_METHOD(RowColumnSum)
		{
			outputLog << "Running UnitTest->RowColumnSum.\n";
			const int N = 256;
			const int iter = 20;
			ArrayXXf X(N, N);
			ArrayXXd Xd(N, N);
			float Xf[N][N];
			float XfRow[N] = {}, XfCol[N] = {};
			ArrayXd resRowDouble(N), resColDouble(N);
			ArrayXf resRow(N), forResRow(N), forResRow2(N), resVec(N), resVecRow(N), resRowOnes(N);
			ArrayXXf resCol(1, N), forResCol(1, N);
			std::vector<ArrayXf> resMat;
			resMat.resize(N);
			for (auto& resmat : resMat) { resmat.resize(N); resmat.setRandom(); }
			float row1 = 0.f;
			float col1 = 0.f;
			resRowDouble.setZero();
			resColDouble.setZero();
			resRow.setZero();
			forResRow.setZero();
			forResRow2.setZero();
			resCol.setZero();
			resVec.setZero();
			resVecRow.setZero();
			forResCol.setZero();
			double durationRow = 1e10, durationCol = 1e10, durationForRow = 1e10, durationForCol = 1e10, durationForRow2 = 1e10, durationRow1 = 1e10, durationCol1 = 1e10, durationVec = 1e10, durationVecRow = 1e10, durationRowMult = 1e10, durationColMult = 1e10;
			double durationRowDouble = 1e10, durationColDouble = 1e10, durationRowFloat = 1e10, durationColFloat = 1e10;
			double durationRowOnes = 1e10;
			for (auto i = 0; i < iter; i++)
			{
				X.setRandom();
				Xd.setRandom();
				for (auto j = 0; j < N; j++) { for (auto k = 0; k < N; k++) { Xf[j][k] = X(j, k); } }

				// rowwise sum
				auto start = std::chrono::steady_clock::now();
				resRow += X.rowwise().sum();
				auto end = std::chrono::steady_clock::now();
				auto time = std::chrono::duration<double, std::micro>(end - start).count();
				durationRow = std::min(durationRow, time);

				// colwise sum
				start = std::chrono::steady_clock::now();
				resCol += X.colwise().sum();
				end = std::chrono::steady_clock::now();
				time = std::chrono::duration<double, std::micro>(end - start).count();
				durationCol = std::min(durationCol, time);

				// for-loop rowwise sum
				start = std::chrono::steady_clock::now();
				for (auto j = 0;j < N;j++) { forResRow(j) += X.row(j).sum(); }
				end = std::chrono::steady_clock::now();
				time = std::chrono::duration<double, std::micro>(end - start).count();
				durationForRow = std::min(durationForRow, time);

				// for-loop rowwise sum 2
				start = std::chrono::steady_clock::now();
				for (auto j = 0;j < N;j++) { forResRow2 += X.col(j); }
				end = std::chrono::steady_clock::now();
				time = std::chrono::duration<double, std::micro>(end - start).count();
				durationForRow2 = std::min(durationForRow2, time);

				// for-loop colwise sum
				start = std::chrono::steady_clock::now();
				for (auto j = 0;j < N;j++) { forResCol(j) += X.col(j).sum(); }
				end = std::chrono::steady_clock::now();
				time = std::chrono::duration<double, std::micro>(end - start).count();
				durationForCol = std::min(durationForCol, time);

				// sum 1-dim row
				start = std::chrono::steady_clock::now();
				for (auto j = 0;j < N;j++) { row1 += resRow.sum(); }
				end = std::chrono::steady_clock::now();
				time = std::chrono::duration<double, std::micro>(end - start).count();
				durationRow1 = std::min(durationRow1, time);

				// sum 1-dim col
				start = std::chrono::steady_clock::now();
				for (auto j = 0;j < N;j++) { col1 += resCol.sum(); }
				end = std::chrono::steady_clock::now();
				time = std::chrono::duration<double, std::micro>(end - start).count();
				durationCol1 = std::min(durationCol1, time);

				// sum 1-dim std:vector
				start = std::chrono::steady_clock::now();
				for (auto j = 0;j < N;j++) { resVec(j) += resMat[j].sum(); }
				end = std::chrono::steady_clock::now();
				time = std::chrono::duration<double, std::micro>(end - start).count();
				durationVec = std::min(durationVec, time);

				// sum 1-dim std:vector over the vectors
				start = std::chrono::steady_clock::now();
				for (auto j = 0;j < N;j++) { resVecRow += resMat[j]; }
				end = std::chrono::steady_clock::now();
				time = std::chrono::duration<double, std::micro>(end - start).count();
				durationVecRow = std::min(durationVecRow, time);

				//multiply a row
				start = std::chrono::steady_clock::now();
				for (auto j = 0;j < N;j++) { X.row(j) *= resCol; }
				end = std::chrono::steady_clock::now();
				time = std::chrono::duration<double, std::micro>(end - start).count();
				durationRowMult = std::min(durationRowMult, time);

				//multiply a column
				start = std::chrono::steady_clock::now();
				for (auto j = 0;j < N;j++) { X.col(j) *= resRow; }
				end = std::chrono::steady_clock::now();
				time = std::chrono::duration<double, std::micro>(end - start).count();
				durationColMult = std::min(durationColMult, time);

				// double rowwise sum
				start = std::chrono::steady_clock::now();
				resRowDouble += Xd.rowwise().sum();
				end = std::chrono::steady_clock::now();
				time = std::chrono::duration<double, std::micro>(end - start).count();
				durationRowDouble = std::min(durationRowDouble, time);

				// double colwise sum
				start = std::chrono::steady_clock::now();
				resColDouble += Xd.colwise().sum();
				end = std::chrono::steady_clock::now();
				time = std::chrono::duration<double, std::micro>(end - start).count();
				durationColDouble = std::min(durationColDouble, time);

				// rowwise sum
				start = std::chrono::steady_clock::now();
				for (auto j = 0; j < N; j++) {
					for (auto k = 0; k < N; k++) { XfRow[j] += Xf[k][j]; }
				}
				end = std::chrono::steady_clock::now();
				time = std::chrono::duration<double, std::micro>(end - start).count();
				durationRowFloat = std::min(durationRowFloat, time);

				// rowwise sum
				start = std::chrono::steady_clock::now();
				for (auto j = 0; j < N; j++) {
					for (auto k = 0; k < N; k++) { XfCol[j] += Xf[j][k]; }
				}
				end = std::chrono::steady_clock::now();
				time = std::chrono::duration<double, std::micro>(end - start).count();
				durationColFloat = std::min(durationColFloat, time);

				// rowwise sum
				start = std::chrono::steady_clock::now();
				resRowOnes.matrix() += X.matrix() * VectorXf::Ones(N);
				end = std::chrono::steady_clock::now();
				time = std::chrono::duration<double, std::micro>(end - start).count();
				durationRowOnes = std::min(durationRowOnes, time);
			}

			for (auto i = 1; i < N; i++) { XfRow[0] += XfRow[i]; }
			outputLog << "Eigen rowwise Time: " << durationRow << " us.\n";
			outputLog << "Eigen columnwise Time: " << durationCol << " us.\n";
			outputLog << "for-loop rowwise Time: " << durationForRow << " us.\n";
			outputLog << "for-loop rowwise2 Time: " << durationForRow2 << " us.\n";
			outputLog << "for-loop colwise Time: " << durationForCol << " us.\n";
			outputLog << "Eigen row 1dim Time: " << durationRow1 << " us.\n";
			outputLog << "Eigen col 1dim Time: " << durationCol1 << " us.\n";
			outputLog << "std::vector rowwise Time: " << durationVecRow << " us.\n";
			outputLog << "std::vector Time: " << durationVec << " us.\n";
			outputLog << "Mult Row Time: " << durationRowMult << " us.\n";
			outputLog << "Mult Col Time: " << durationColMult << " us.\n";
			outputLog << "Eigen double precision rowwise Time: " << durationRowDouble << " us.\n";
			outputLog << "Eigen double precision columnwise Time: " << durationColDouble << " us.\n";
			outputLog << "Float array rowwise Time: " << durationRowFloat << " us.\n";
			outputLog << "Float array columnwise Time: " << durationColFloat << " us.\n";
			outputLog << "Eigen rowwise mult with ones Time: " << durationRowOnes << " us.\n";
			outputLog << "output error is : " << resRow.sum() << ", " << forResRow.sum() << ", " << XfRow[0] << std::endl;
			Assert::IsTrue(resRow.abs2().sum() == forResRow.abs2().sum());
			Assert::IsTrue(resCol.abs2().sum() == forResCol.abs2().sum());
			float error = fabs(row1 - col1) / fabs(row1); // strangely enough these are not exactly the same...
			Assert::IsTrue(error*error < 1e-10f);

		}
#endif // #if NDEBUG
	};
}

namespace Algorithms
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
		bool flag = algoStreaming.Initialize(static_cast<int>(xTime.rows() / 2) , static_cast<int>(xTime.cols()), sampleRate);
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
		flag = algoStreaming.Initialize(static_cast<int>(xTime.rows() * 1.25f), static_cast<int>(xTime.cols() ), sampleRate);
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

	TEST_CLASS(BeamformerAdaptiveTest)
	{
		class BeamformerAdaptiveDebug : public BeamformerAdaptive
		{
		public:
			using BeamformerAdaptive::GetPrivateData;
			using BeamformerAdaptive::CalculateFilter;
		};
		TEST_METHOD(Interface)
		{
			outputLog << "Running BeamformerAdaptiveTest->Interface.\n";
			BeamformerAdaptive beamformer;
			auto c = beamformer.GetCoefficients();
			ArrayXXcf input1(c.NBands, c.NChannels);
			bool input2 = true;
			input1.setRandom();
			I::BeamformerAdaptive input = { input1, input2 };
			ArrayXcf output(c.NBands);
			Assert::IsTrue(AlgorithmInterfaceTest<BeamformerAdaptive>(input, output));
		}
#if NDEBUG
		TEST_METHOD(CovarianceMatrices)
		{
			outputLog << "Running BeamformerAdaptiveTest->CovarianceMatrices: Testing different covariance scenarios.\n";

			BeamformerAdaptiveDebug beamformer;
			auto c = beamformer.GetCoefficients();
			c.FilterbankRate = 125.f;
			c.NBands = 257;
			beamformer.Initialize(c);
			auto p = beamformer.GetParameters();
			p.FilterUpdateTConstant = 0.25f; // set time constant of covariance matrix updates to .25 second
			beamformer.SetParameters(p);

			// run for 3 seconds with random noise
			int Nend = int(3.f * c.FilterbankRate);
			float NoiseLevelL = 0.6f;
			float NoiseLevelR = 0.2f;
			ArrayXcf y = ArrayXcf::Zero(c.NBands);
			ArrayXXcf noise;
			ArrayXcf noiseR;
			for (auto n = 0; n < Nend; n++)
			{
				noise = ArrayXXcf::Random(c.NBands, c.NChannels);
				noise.col(0) *= NoiseLevelL;
				noise.col(1) *= NoiseLevelR;
				beamformer.Process({ noise, false }, y);
			}

			auto Data = beamformer.GetPrivateData();
			// Run again for 3 seconds and measure error
			float rnL = 0.f;
			float rnR = 0.f;
			for (auto n = 0; n < Nend; n++)
			{
				noise = ArrayXXcf::Random(c.NBands, c.NChannels);
				noise.col(0) *= NoiseLevelL;
				noise.col(1) *= NoiseLevelR;
				beamformer.Process({ noise, false }, y);
				for (auto i = 0;i < c.NBands;i++)
				{
					rnL += std::abs(Data->Rn[i](0, 0));
					rnR += std::abs(Data->Rn[i](1, 1));
				}
			}
			rnL /= (c.NBands * Nend);
			rnR /= (c.NBands * Nend);
			float errorL = abs(rnL - .666f * NoiseLevelL * NoiseLevelL) / (.666f * NoiseLevelL * NoiseLevelL);
			float errorR = abs(rnR - .666f * NoiseLevelR * NoiseLevelR) / (.666f * NoiseLevelR * NoiseLevelR);

			outputLog << "Left mean error is: " << errorL << ".\n";
			outputLog << "Right mean error is: " << errorR << ".\n";
			Assert::IsTrue(errorL < .05f);
			Assert::IsTrue(errorR < .05f);

			// run again, and every 3rd sample add a large correlated transient
			float rnLC = 0.f;
			float rnRC = 0.f;
			for (auto n = 0; n < Nend; n++)
			{

				if (n % 3 == 0) {
					noise.col(0) = 10 * NoiseLevelL * ArrayXcf::Random(c.NBands);
					noise.col(1) = noise.col(0);
					beamformer.Process({ noise, true }, y);
				}
				else {
					noise.col(0) = NoiseLevelL * ArrayXcf::Random(c.NBands);
					noise.col(1) = NoiseLevelR * ArrayXcf::Random(c.NBands);
					beamformer.Process({ noise, false }, y);
				}

				for (auto i = 0;i < c.NBands;i++)
				{
					rnLC += std::abs(Data->Rn[i](0, 0));
					rnRC += std::abs(Data->Rn[i](1, 1));
				}
			}
			rnLC /= (c.NBands * Nend);
			rnRC /= (c.NBands * Nend);
			errorL = abs(rnLC - rnL) / rnL;
			errorR = abs(rnRC - rnR) / rnR;

			outputLog << "Left mean error is: " << errorL << ".\n";
			outputLog << "Right mean error is: " << errorR << ".\n";
			Assert::IsTrue(errorL < .05f);
			Assert::IsTrue(errorR < .05f);
			ArrayXf CrossX = ArrayXf::Zero(c.NBands);
			ArrayXf CrossN = ArrayXf::Zero(c.NBands);
			for (auto i = 0;i < c.NBands;i++)
			{
				CrossX(i) = std::abs(Data->Rx[i](1, 0)) / std::sqrt(std::abs(Data->Rx[i](0, 0)) * std::abs(Data->Rx[i](1, 1)));
				CrossN(i) = std::abs(Data->Rn[i](1, 0)) / std::sqrt(std::abs(Data->Rn[i](0, 0)) * std::abs(Data->Rn[i](1, 1)));
				Assert::IsTrue((CrossX(i) > 0.9f));
				Assert::IsTrue((CrossN(i) < 0.4f));
			}
			for (auto i = 0;i < c.NBands;i++)
			{
				Assert::IsTrue(std::abs(Data->Rx[i](1, 0)) > 10 * std::abs(Data->Rn[i](1,0)));
			}
			outputLog << "Test succesful" << std::endl;
		}
#endif //#if NDEBUG
		TEST_METHOD(Disabled)
		{
			using namespace std::literals;

			outputLog << "Running BeamformerAdaptiveTest->Disabled.\n";
			BeamformerAdaptive beamformer;
			auto c = beamformer.GetCoefficients();
			beamformer.Disable();

			const auto numberOfFrames = 10;
			ArrayXcf output = ArrayXcf::Zero(c.NBands);
			for (int i = 0; i < numberOfFrames; i++)
			{
				ArrayXXcf input = ArrayXXcf::Random(c.NBands, c.NChannels) + 1.0if * ArrayXXcf::Random(c.NBands, c.NChannels);
				beamformer.Process({ input,true }, output);
				Assert::IsTrue((output - input.col(0)).abs().sum() == 0.0);
			}
		}

		TEST_METHOD(FilterCalculation)
		{
			using namespace std::literals;

			outputLog << "Running BeamformerAdaptiveTest->FilterCalculation: Testing the Filter calculation.\n";
			BeamformerAdaptiveDebug beamformer;
			auto c = beamformer.GetCoefficients();
			beamformer.Initialize(c);

			auto Data = beamformer.GetPrivateData();
			for (auto i = 0;i < c.NBands;i++)
			{
				Data->Rx[i](0, 0) = .5f;
				Data->Rx[i](1, 1) = .5f;
				Data->Rx[i](1, 0) = .0f;

				Data->Rn[i](0, 0) = .2f;
				Data->Rn[i](1, 1) = .2f;
				Data->Rn[i](1, 0) = .0f;

				Data->Rx[i] -= Data->Rn[i];
			}

			for (int i = 0; i < c.NBands; i++)
			{
				Data->CurrentBand = i;
				beamformer.CalculateFilter();
			}

			Assert::IsTrue(Data->Filter.abs2().sum() < 1e-10);


			for (auto i = 0;i < c.NBands;i++)
			{
				Data->Rx[i](0, 0) = .7f;
				Data->Rx[i](1, 1) = .7f;
				Data->Rx[i](1, 0) = .7f;

				Data->Rn[i](0, 0) = .5f;
				Data->Rn[i](1, 1) = .5f;
				Data->Rn[i](1, 0) = .5f;

				Data->Rx[i] -= Data->Rn[i];
			}

			for (int i = 0; i < c.NBands; i++)
			{
				Data->CurrentBand = i;
				beamformer.CalculateFilter();
			}
			auto err = Data->Filter.abs2().sum();
			//Assert::IsTrue(beamformer.D.Filter.abs2().sum() < 1e-10);

			Assert::IsTrue((Data->Filter.col(0) - 1.f).abs2().sum() < 1e-10);
			Assert::IsTrue(Data->Filter.col(1).abs2().sum() < 1e-10);

			for (auto i = 0;i < c.NBands;i++)
			{
				Data->Rx[i](0, 0) = .7792f;
				Data->Rx[i](1, 1) = .934f;
				Data->Rx[i](1, 0) = -.8213 + .2307i;

				Data->Rn[i](0, 0) = .4694f;
				Data->Rn[i](1, 1) = .0119f;
				Data->Rn[i](1, 0) = -.0178 - .0726i;

			}

			for (int i = 0; i < c.NBands; i++)
			{
				Data->CurrentBand = i;
				beamformer.CalculateFilter();
			}

			float error = (Data->Filter.col(0) - (.0253 - .141i)).abs2().mean();
			Assert::IsTrue(error < 1e-6);
			error = (Data->Filter.col(1) - (-.8222 - .3647i)).abs2().mean();
			Assert::IsTrue(error < 1e-6);

			outputLog << "Test succesful" << std::endl;
		}

		TEST_METHOD(FreezeOrAutomatic)
		{
			using namespace std::literals;

			outputLog << "Running BeamformerAdaptiveTest->FreezeOrAutomatic.\n";
			BeamformerAdaptive beamformer;
			auto c = beamformer.GetCoefficients();
			auto p = beamformer.GetParameters();
			p.SpeechDecision = p.FreezeUpdate;
			beamformer.Initialize(c);
			beamformer.SetParameters(p);

			const auto numberOfFrames = 100;
			ArrayXcf output = ArrayXcf::Zero(c.NBands);
			for (int i = 0; i < numberOfFrames; i++)
			{
				ArrayXXcf input = ArrayXXcf::Random(c.NBands, c.NChannels) + 1.0if * ArrayXXcf::Random(c.NBands, c.NChannels);
				beamformer.Process({ input, true }, output);
				Assert::IsTrue((output - input.col(0)).abs().sum() == 0.0);
			}

			p.SpeechDecision = p.Automatic;
			beamformer.SetParameters(p);
			ArrayXcf outputAutomatic = ArrayXcf::Zero(c.NBands);
			for (int i = 0; i < 10; i++)
			{
				ArrayXXcf input = ArrayXXcf::Random(c.NBands, c.NChannels) + 1.0if * ArrayXXcf::Random(c.NBands, c.NChannels);
				beamformer.Process({ input, true }, outputAutomatic);
			}
			Assert::IsTrue((output - outputAutomatic).abs().sum() != 0.0);
		}

#if NDEBUG
		TEST_METHOD(WienerGain)
		{
			outputLog << "Running BeamformerAdaptiveTest->WienerGain.\n";
			srand(1);

			BeamformerAdaptive beamformer;
			VoiceActivationDetection VAD;
			auto cBF = beamformer.GetCoefficients();
			auto cVAD = VAD.GetCoefficients();
			cVAD.FilterbankRate = cBF.FilterbankRate;
			cVAD.NBands = cBF.NBands;
			cVAD.NChannels = cBF.NChannels;
			beamformer.Initialize(cBF);
			VAD.Initialize(cVAD);

			

			const int noFrames = static_cast<int>(30 * cBF.FilterbankRate); // 30 seconds
			ArrayXXcf input1(cBF.NBands, cBF.NChannels);
			ArrayXcf output(cBF.NBands);
			bool voiceActivity;

			for (auto i = 0;i < noFrames;i++)
			{
				input1 = ArrayXXcf::Random(cBF.NBands, cBF.NChannels);
				VAD.Process(input1, voiceActivity);
				beamformer.Process({ input1, voiceActivity }, output);

			}

			// measure output with SpeechDistortionRatio = 1
			auto pBF = beamformer.GetParameters();
			pBF.SpeechDistortionRatio = 1.f;
			beamformer.SetParameters(pBF);
			
			srand(1);
			ArrayXf powerOut1(cBF.NBands);
			powerOut1.setZero();
			for (auto i = 0;i < noFrames;i++)
			{
				input1 = ArrayXXcf::Random(cBF.NBands, cBF.NChannels);
				VAD.Process(input1, voiceActivity);
				beamformer.Process({ input1, voiceActivity }, output);
				powerOut1 += output.abs2() / noFrames;
			}

			// measure output with SpeechDistortionRatio = 0
			pBF.SpeechDistortionRatio = 0.f;
			beamformer.SetParameters(pBF);

			srand(1);
			ArrayXf powerOut0(cBF.NBands);
			powerOut0.setZero();
			for (auto i = 0;i < noFrames;i++)
			{
				input1 = ArrayXXcf::Random(cBF.NBands, cBF.NChannels);
				VAD.Process(input1, voiceActivity);
				beamformer.Process({ input1, voiceActivity }, output);
				powerOut0 += output.abs2() / noFrames;
			}

			float attenuationdB = 10 * log10f((powerOut0 / powerOut1).mean());
			outputLog << "Attenuation: " << attenuationdB << " dB\n";
			outputLog << "Test successful." << std::endl;
		}
#endif // #if NDEBUG
		TEST_METHOD(ProcessComplexArrays)
		{
			outputLog << "Running BeamformerAdaptiveTest->ProcessComplexArrays: Testing the input/output.\n";
			srand(1); // reset random generator so calls to random generator are deterministic

			BeamformerAdaptive beamformer;
			VoiceActivationDetection voiceActivationDetection;
			auto c = beamformer.GetCoefficients();
			c.NBands = 129;
			c.FilterbankRate = 125;
			c.NChannels = 2;
			beamformer.Initialize(c);
			auto cVAD = voiceActivationDetection.GetCoefficients();
			cVAD.FilterbankRate = 125;
			cVAD.NBands = 129;
			cVAD.NChannels = 2;
			voiceActivationDetection.Initialize(cVAD);

			AudioFile<float> audioFile;
			audioFile.setNumChannels(2);
			audioFile.setBitDepth(24);
			audioFile.setSampleRate(16000);

			const int noFrames = 500;

			ArrayXcf outputx = ArrayXcf::Zero(c.NBands);
			ArrayXcf outputArray = ArrayXcf::Zero(c.NBands * noFrames);


			outputLog << "Running Beamformer....\n";
			bool voiceActivity;
			for (auto i = 0; i < noFrames; i++)
			{
				ArrayXXcf inputx = ArrayXXcf::Random(c.NBands, c.NChannels);
				voiceActivationDetection.Process(inputx, voiceActivity);
				beamformer.Process({ inputx , voiceActivity }, outputx);

				for (auto j = 0; j < c.NBands; j++)
				{
					audioFile.samples[0].push_back(outputx(j).real());
					audioFile.samples[1].push_back(outputx(j).imag());
				}
			}

			outputLog << "Saving file.\n";
			audioFile.save("../../../../Temp/BeamformerAdaptive.wav", AudioFileFormat::Wave); // truncate to 24bits 
			Assert::IsTrue(audioFile.load("../../../../Temp/BeamformerAdaptive.wav")); // load again

			AudioFile<float> audioFileRef;
			audioFileRef.setNumChannels(2);
			audioFileRef.setBitDepth(24);
			audioFileRef.setSampleRate(16000);
			Assert::IsTrue(audioFileRef.load("../../../../Datasets/wav/BeamformerAdaptive.wav")); // compare with reference

			outputLog << "Calculating error...\n";
			auto N = audioFile.getNumSamplesPerChannel();
			float Error = 0.f;
			for (auto i = 0; i < N; i++)
			{
				Error += abs(audioFileRef.samples[0][i] - audioFile.samples[0][i]) + abs(audioFileRef.samples[1][i] - audioFile.samples[1][i]);
			}
			Error /= (c.NBands * noFrames * 2);
			outputLog << "Error: " << Error << std::endl;
			Assert::IsTrue(Error < 1e-10);

			outputLog << "Test successful." << std::endl;
		}
	};

	TEST_CLASS(CubicSplineTest)
	{
		TEST_METHOD(Interface)
		{
			outputLog << "Running CubicSplineTest->Interface.\n";
			CubicSpline spline;
			int N = 11; // Number of known points
			int Nsplines = 2; // Number of splines to interpolate in one function call
			int Ndesired = 1001; // desired number of points to estimate
			ArrayXXf input1(N, Nsplines);
			ArrayXXf output1(N, Nsplines);
			ArrayXXf xDesired(Ndesired, Nsplines);
			input1 = VectorXf::LinSpaced(N, 0, 1000).array().replicate(1,Nsplines); // known x values
			output1.setRandom(); // random known y values
			xDesired = VectorXf::LinSpaced(Ndesired, 0, 1000).array().replicate(1,Nsplines); // desired x values
			I::CubicSpline input = { input1, output1, xDesired };

			ArrayXXf yDesired(Ndesired, Nsplines);
			Assert::IsTrue(AlgorithmInterfaceTest<CubicSpline>(input, yDesired));
		}

		TEST_METHOD(CheckCalculation)
		{
			CubicSpline spline;
			spline.Initialize();
			int N = 3;
			int NDesired = 10;
			ArrayXf input1(N), output1(N), xDesired(NDesired), yDesired(NDesired), yRef(NDesired);
			input1 << 1.f, 2.f, 3.f;
			output1 << 2.f, 6.f, 4.3f;
			xDesired << 1.f, 1.2f, 1.67f, 1.8f, 2.1f, 2.45f, 2.6f, 2.7f, 2.9f, 3.f;
			I::CubicSpline input = { input1, output1, xDesired };
			spline.Process(input, yDesired);

			yRef << 2.f, 2.8912f, 5.1022f, 5.5648f, 6.0609f, 5.6230f, 5.2536f, 4.9895f, 4.4956f, 4.3f;
			float error = (yRef - yDesired).abs2().sum() / yRef.abs2().sum();
			outputLog << "Error: " << error << std::endl;
			Assert::IsTrue(error < 1e-5);
		}

		TEST_METHOD(InputReverseOrder)
		{
			CubicSpline spline;
			spline.Initialize();
			int N = 3;
			int NDesired = 10;
			ArrayXf input1(N), output1(N), xDesired(NDesired), yDesired(NDesired), yRef(NDesired);
			input1 << 3.f, 2.f, 1.f;
			output1 << 4.3f, 6.f, 2.f;
			xDesired << 1.f, 1.2f, 1.67f, 1.8f, 2.1f, 2.45f, 2.6f, 2.9f, 2.7f, 3.f;
			I::CubicSpline input = { input1, output1, xDesired };
			spline.Process(input, yDesired);

			yRef << 2.f, 2.8912f, 5.1022f, 5.5648f, 6.0609f, 5.6230f, 5.2536f, 4.4956f, 4.9895f, 4.3f;
			float error = (yRef - yDesired).abs2().sum() / yRef.abs2().sum();
			outputLog << "Error: " << error << std::endl;
			Assert::IsTrue(error < 1e-5);
		}
	};

	TEST_CLASS(DesignIIRMinPhaseTest)
	{
		TEST_METHOD(Interface)
		{
			using namespace std::literals;
			outputLog << "Running DesignIIRMinPhaseTest->Interface.\n";
			DesignIIRMinPhase filterDesigner;
			auto c = filterDesigner.GetCoefficients();
			ArrayXf input(c.NBands);
			input.setRandom();
			ArrayXXf SOS(filterDesigner.GetNSOS(c.N), 6);
			float gain;
			O::DesignIIRMinPhase output{ SOS, gain };
			Assert::IsTrue(AlgorithmInterfaceTest<DesignIIRMinPhase>(input, output));
		}

		TEST_METHOD(CheckCalculation)
		{
			using namespace std::literals;
			outputLog << "Running DesignIIRMinPhaseTest->CheckCalculation.\n";
			DesignIIRMinPhase filterDesigner;
			auto c = filterDesigner.GetCoefficients();
			c.N = 5;
			c.NBands = 17;
			c.WeightType = c.Linear;
			filterDesigner.Initialize(c);
			ArrayXf input(c.NBands);
			input << 1.4560f, 2.2489f, 0.9559f, 5.0838f, 1.2221f, 2.3314f, 1.9187f, 1.5807f, 1.9052f, 1.8399f, 6.1568f, 1.3541f, 0.5104f, 1.9502f, 1.5755f, 4.1232f, 5.1997f;
			ArrayXXf SOS(filterDesigner.GetNSOS(c.N), 6);
			float gain;
			O::DesignIIRMinPhase output{ SOS, gain };
			filterDesigner.Process(input, output);
			ArrayXXf SOSref(filterDesigner.GetNSOS(c.N), 6);
			SOSref << 1.f, 1.04831f, 0.842527f, 1.f, -1.50516f, 0.884276f, 1.f, -1.34835f, 0.751896f, 1.f, 0.759416f, 0.870153f, 1.f, 0.345899f, 0.f, 1.f, 0.850303f, 0.f;
			float errorSOS = (SOS - SOSref).abs2().sum() / SOSref.abs2().sum();
			outputLog << "SOS relative error: " << errorSOS << "\n";
			float errorGain = std::fabs(gain - 1.89722f);
			errorGain = errorGain * errorGain / (1.89722f*1.89722f);
			outputLog << "gain relative error: " << errorGain << std::endl;
			outputLog << "SOS: " << SOS << std::endl;
			Assert::IsTrue(errorSOS < 1e-10f);
			Assert::IsTrue(errorGain < 1e-10f);
		}
	};

	TEST_CLASS(DesignIIRNonParametricTest)
	{
		TEST_METHOD(Interface)
		{
			using namespace std::literals;
			outputLog << "Running DesignIIRNonParametricTest->Interface.\n";
			DesignIIRNonParametric filterDesigner;
			auto c = filterDesigner.GetCoefficients();
			ArrayXf frequencies(5);
			ArrayXf gaindB(5);
			frequencies << 80, 500, 1000, 2400, 4000;
			gaindB << 5, 10, -10, 0, 7;
			I::DesignIIRNonParametric input = { frequencies, gaindB };
			ArrayXXf SOS(filterDesigner.FilterDesigner.GetNSOS(c.N), 6);
			float gain;
			O::DesignIIRNonParametric output{ SOS, gain };
			Assert::IsTrue(AlgorithmInterfaceTest<DesignIIRNonParametric>(input, output));
		}
	};

	TEST_CLASS(NoiseEstimationTest)
	{
	public:
		TEST_METHOD(Interface)
		{
			outputLog << "Running NoiseEstimationTest->Interface.\n";
			NoiseEstimationSPP noiseEstimation;
			auto c = noiseEstimation.GetCoefficients();
			ArrayXXf input(c.NBands, c.NChannels);
			input.setRandom();
			ArrayXXf output1(c.NBands, c.NChannels), output2(c.NBands, c.NChannels);
			O::NoiseEstimationSPP output{ output1, output2 };
			Assert::IsTrue(AlgorithmInterfaceTest<NoiseEstimationSPP>(input, output));
		}

		TEST_METHOD(NoiseStationaryConvergence)
		{
			outputLog << "Running NoiseEstimationTest->NoiseStationaryConvergence.\n";

			srand(1);  // reset random generator to get consistent results
			NoiseEstimationSPP noiseEstimation;
			auto c = noiseEstimation.GetCoefficients();
			c.NChannels = 1;
			noiseEstimation.Initialize(c);
			
			ArrayXf powerInput = (ArrayXf::Random(c.NBands) + 1.f) * .5f; // array of random numbers between 0 and 1
			ArrayXf activity = ArrayXf::Zero(c.NBands);
			ArrayXf powerNoise(c.NBands);
			powerNoise.setZero();

			// run once and check that Activity == 1
			noiseEstimation.Process(powerInput, { powerNoise, activity });
			float ADError = 1.f / c.NBands * ((activity - 1.f).abs()).sum();
			outputLog << "Activity detection mean error is: " << ADError << ".\n";
			Assert::IsTrue(ADError < 1e-10);

			// run algorithm for 3 seconds
			int Nend = int(3.f * c.FilterbankRate); // number of samples corresponding to 3 seconds
			for (auto n = 0; n < Nend; n++) {
				noiseEstimation.Process(powerInput, { powerNoise, activity });
			}

			// calculate error
			float Error = 1.f / c.NBands * ((powerNoise - powerInput).abs2()).sum();
			outputLog << "Noise mean-square error is: " << Error << ".\n";
			Assert::IsTrue(Error < 1e-10);

			// check Activity < .1
			ADError = 1.f / c.NBands * activity.sum();
			outputLog << "Activity mean error is: " << ADError << ".\n";
			Assert::IsTrue(ADError < .1f);

			outputLog << "Test succesful." << std::endl;
		}
		TEST_METHOD(TestMethod1)
		{
			NoiseEstimationSPP noiseEstimation;
			auto c = noiseEstimation.GetCoefficients();
			outputLog << c.FilterbankRate << std::endl;
			noiseEstimation.Initialize(c);
			auto p = noiseEstimation.GetParameters();
			outputLog << p.ActivityMeanTConstant << std::endl;
			noiseEstimation.SetParameters(p);
			auto m = noiseEstimation.GetAllocatedMemorySize();

			auto t = noiseEstimation.GetParametersAll();
			//t.Parameters.ActivityMeanTConstant = 1.f;
			//noiseEstimation.SetParametersAll(t);
			//auto coef = noiseEstimation.GetCoefficientsAll();
			//auto s = noiseEstimation.GetSetupAll();
			//noiseEstimation.SetSetupAll(s);
			//outputLog << sizeof(s) << std::endl;
			////outputLog << t << std::endl;
			////t.P = 1;
			//VectorAlgo<NoiseEstimationSPP> test;
			//
			//test.resize(10);
			//test.Initialize();
			//outputLog << test[3].GetParameters().ActivityMeanTConstant << std::endl;
			//test[3].GetParameters().ActivityMeanTConstant = 10;
			//outputLog << test[3].GetParameters().ActivityMeanTConstant << std::endl;

			noiseEstimation.Reset();
			noiseEstimation.Initialize();
			outputLog << noiseEstimation.GetAllocatedMemorySize() << std::endl;
			//noiseEstimation.fft.resize(10);

		}

	};


	TEST_CLASS(GainCalculationTest)
	{
		TEST_METHOD(Interface)
		{
			outputLog << "Running GainCalculationTest->Interface.\n";
			GainCalculationSimple gainCalculation;
			auto c = gainCalculation.GetCoefficients();
			ArrayXXf input(100,10);
			input.setRandom();
			ArrayXXf output(100,10);
			auto flag = AlgorithmInterfaceTest<GainCalculationSimple>(input, output);
			Assert::IsTrue(flag);
		}
	};

	TEST_CLASS(InterpolationTemporalTest)
	{
		TEST_METHOD(Interface)
		{
			InterpolationTemporal timeStretch;
			auto c = timeStretch.GetCoefficients();
			ArrayX2f xEnergy(c.NBands,2);
			ArrayX2f xPhase(c.NBands,2);
			float point = 0.2f;
			xEnergy.setRandom().abs();
			xPhase.setRandom();
			I::InterpolationTemporal input = { xEnergy, xPhase, point };
			ArrayXcf output(c.NBands);
			auto flag = AlgorithmInterfaceTest<InterpolationTemporal>(input, output);
			Assert::IsTrue(flag);

		}
	};

	TEST_CLASS(MinPhaseSpectrumTest)
	{
		TEST_METHOD(Interface)
		{
			outputLog << "Running MinPhaseSpectrumTest->Interface.\n";
			MinPhaseSpectrum mps;
			auto c = mps.GetCoefficients();
			auto p = mps.GetParameters();
			auto nChannels = 2;
			ArrayXXf input(c.NBands, nChannels);
			input.setRandom();
			ArrayXXcf output(c.NBands, nChannels);
			auto flag = AlgorithmInterfaceTest<MinPhaseSpectrum>(input, output);
			Assert::IsTrue(flag);
		}

		TEST_METHOD(CheckCalculation)
		{
			using namespace std::complex_literals;

			MinPhaseSpectrum mps;
			int N = 17;
			auto c = mps.GetCoefficients();
			c.NBands = N;			
			mps.Initialize(c);
			ArrayXf input = ArrayXf::LinSpaced(N, 1.f, static_cast<float>(N));
			ArrayXcf output(N);
			mps.Process(input, output);
			ArrayXcf outRef(N);
			outRef << 1.f, 1.2739f + 1.5418if, 2.0368f + 2.2026if, 2.2561f + 3.3030if, 3.1422f + 3.8893if, 3.5936f + 4.8048if, 4.6479f + 5.2342if, 5.3500f + 5.9479if, 6.5712f + 6.1498if, 7.5192f + 6.5926if, 8.8878f + 6.4812if, 10.0579f + 6.5451if, 11.5362f + 5.9931if, 12.8699f + 5.5105if, 14.3793f + 4.2704if, 15.7290f + 2.9322if, 17.0000f;
			float error = (output - outRef).abs2().sum() / outRef.abs2().sum();
			outputLog << "Error: " << error << std::endl;
			Assert::IsTrue(error < 1e-10f);
		}
	};

	TEST_CLASS(InterpolationCubicTest)
	{
		TEST_METHOD(Interface)
		{
			ArrayXf inputp(4);
			inputp << 1.f, 3.2f, 2.6f, 2.9f;
			ArrayXf d(9);
			d << 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f;
			I::InterpolationCubic input = { inputp, d };
			ArrayXf output(9);
			auto flag = AlgorithmInterfaceTest<InterpolationCubic>(input, output);
			Assert::IsTrue(flag);
		}

		TEST_METHOD(ValuesTest)
		{
			InterpolationCubic interpolate;
			interpolate.Initialize();
			ArrayXf inputp(4);
			inputp << 1.f, 3.2f, 2.6f, 2.9f;
			ArrayXf d(9);
			d << 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f;
			I::InterpolationCubic input = { inputp, d };
			ArrayXf output(9), ref(9);
			interpolate.Process(input, output);
			ref << 3.24935f, 3.2448f, 3.19745f, 3.1184f, 3.01875f, 2.9096f, 2.80205f, 2.7072f, 2.63615f;
			auto error = (output - ref).abs2().sum();
			outputLog << error << std::endl;
			Assert::IsTrue(error < 1e-10f);
		}

		TEST_METHOD(InterfaceConstantDelay)
		{
			ArrayXXf input(4, 1000);
			input.setRandom();
			ArrayXf output(1000);
			auto flag = AlgorithmInterfaceTest<InterpolationCubicConstantDelay>(input, output);
			Assert::IsTrue(flag);
		}

		TEST_METHOD(ValuesConstantDelayTest)
		{
			InterpolationCubicConstantDelay interpolate;
			interpolate.Initialize();
			auto p = interpolate.GetParameters();
			p.FractionalDelay = 0.2f;
			interpolate.SetParameters(p);
			Array4f inputp;
			inputp << 1.f, 3.2f, 2.6f, 2.9f;
			float output, ref;
			ref = 3.2448f;
			interpolate.Process(inputp, Map<ArrayXf>(&output,1));
			auto error = std::abs(ref - output);
			Assert::IsTrue(error < 1e-10f);

		}
	};

	TEST_CLASS(DesignFIRMinPhaseTest)
	{
		TEST_METHOD(Interface)
		{
			outputLog << "Running DesignFIRMinPhaseTest->Interface.\n";
			DesignFIRMinPhase mpsFIRCalculation;
			auto p = mpsFIRCalculation.GetParameters();
			auto c = mpsFIRCalculation.GetCoefficients();
			auto nBands = c.FilterSize / 2 + 1;
			ArrayXXf input(nBands,2);
			input.setRandom().abs();
			ArrayXXf output(c.FilterSize, 2);
			auto flag = AlgorithmInterfaceTest<DesignFIRMinPhase>(input, output);
			Assert::IsTrue(flag);
		}
	};

	TEST_CLASS(DesignFIRNonParametricTest)
	{
		TEST_METHOD(Interface)
		{
			outputLog << "Running DesignFIRNonParametricTest->Interface.\n";
			DesignFIRNonParametric mpsFIRCalculation;
			auto c = mpsFIRCalculation.GetCoefficients();
			ArrayXf frequencies(5);
			ArrayXf gaindB(5);
			frequencies << 80, 500, 1000, 2400, 4000;
			gaindB << 5, 10, -10, 0, 7;
			I::DesignFIRNonParametric input = { frequencies, gaindB };
			ArrayXf output(c.FilterSize);
			auto flag = AlgorithmInterfaceTest<DesignFIRNonParametric>(input, output);
			Assert::IsTrue(flag);
		}

		TEST_METHOD(CheckCalculation)
		{
			DesignFIRNonParametric filter;
			auto c = filter.GetCoefficients();
			auto p = filter.GetParameters();
			c.FilterSize = 32;
			int N = 4;
			filter.Initialize(c);
			ArrayXf frequencies(N), gaindB(N);
			frequencies << 100.f, 1000.f, 3000.f, 6000.f;
			gaindB << 10.f, -20.f, 5.f, 10.f;
			I::DesignFIRNonParametric input = { frequencies, gaindB };
			ArrayXf output(c.FilterSize);
			filter.Process(input, output);
			ArrayXf outRef(c.FilterSize);
			outRef << 1.47128f, -1.78612f, 0.538379f, 0.580462f, 0.221934f, -0.0496397f, 0.0389704f, 0.161114f, 0.179561f, 0.1761f, 0.151415f, 0.12479f, 0.125126f, 0.135422f, 0.14529f, 0.14383f, 0.115281f, 0.111329f, 0.124507f, 0.112827f, 0.102322f, 0.0966862f, 0.0928112f, 0.0874154f, 0.0812955f, 0.0752866f, 0.069461f, 0.0641043f, 0.059177f, 0.0544371f, 0.0497704f, 0.0452115f;
			float error = (output - outRef).abs2().sum() / outRef.abs2().sum();
			outputLog << "Error: " << error << std::endl;
			Assert::IsTrue(error < 1e-8f);
		}
	};

	TEST_CLASS(NoiseSuppressionTest)
	{
		TEST_METHOD(Interface)
		{
			outputLog << "Running NoiseSuppressionTest->Interface.\n";
			NoiseSuppression noiseSuppression;
			auto c = noiseSuppression.GetCoefficients();
			ArrayXXcf input(c.NBands, c.NChannels);
			input.setRandom();
			ArrayXXcf output(c.NBands, c.NChannels);
			Assert::IsTrue(AlgorithmInterfaceTest<NoiseSuppression>(input, output));
		}
	};

	TEST_CLASS(NonparametricEqualizerTest)
	{
		TEST_METHOD(InterfaceStreaming)
		{
			outputLog << "Running NonparametricEqualizerTest->InterfaceStreaming.\n";
			NonparametricEqualizer equalizer;
			auto c = equalizer.GetCoefficients();
			auto p = equalizer.GetParameters();
			ArrayXXf input(c.BufferSize, c.NChannels);
			input.setRandom();
			ArrayXf frequencies(5);
			ArrayXf gaindB(5);
			frequencies << 80, 500, 1000, 2400, 4000;
			gaindB << 5, 10, -10, 0, 7;
			I::NonparametricEqualizerPersistent inputX = { frequencies, gaindB};
			ArrayXXf output(c.BufferSize, c.NChannels);
			float sampleRate = c.SampleRate;
			auto flag = AlgorithmInterfaceStreamingTest<NonparametricEqualizerStreaming>(input, output, sampleRate, inputX);
			Assert::IsTrue(flag);
		}
	};

	TEST_CLASS(PipelinePreProcessingTest)
	{
		TEST_METHOD(InterfaceStreaming)
		{
			outputLog << "Running PipelinePreProcessingTest->Interface.\n";
			PipelinePreProcessing PPP;
			auto c = PPP.GetCoefficients();
			ArrayXXf input(c.BufferSize, c.NChannels);
			input.setRandom();
			ArrayXf output(c.BufferSize);
			auto flag = AlgorithmInterfaceStreamingTest<PipelinePreProcessingStreaming>(input, output, c.SampleRate);
			Assert::IsTrue(flag);
		}
	};

	TEST_CLASS(PitchShiftTest)
	{
		TEST_METHOD(InterfaceStreaming)
		{
			PitchShift timeStretch;
			auto c = timeStretch.GetCoefficients();
			ArrayXXf input(c.BufferSize, 1);
			input.setRandom();
			ArrayXXf output(c.BufferSize, 1);
			float sampleRate = 48e3f;
			auto flag = AlgorithmInterfaceStreamingTest<PitchShiftStreaming>(input, output, sampleRate);
			Assert::IsTrue(flag);
		}

	};

	TEST_CLASS(TonalDetectionTest)
	{
		TEST_METHOD(Interface)
		{
			TonalDetection tonalDetector;
			auto c = tonalDetector.GetCoefficients();
			ArrayXf input(c.NBands);
			input.setZero();
			ArrayXf output(c.NBands);
			auto flag = AlgorithmInterfaceTest<TonalDetection>(input, output);
			Assert::IsTrue(flag);
		}

		TEST_METHOD(CheckCalculation)
		{
			TonalDetection tonalDetector;
			auto setup = tonalDetector.GetSetup();
			setup.Coefficients.FilterbankRate = 5512.5f;
			setup.Coefficients.FilterRelativeRange = 0.1875f;
			setup.Coefficients.NBands = 17;
			setup.Parameters.TonalTC = 0.02f;
			setup.Parameters.TonalTransientRatio = 0.262853538008622f;
			tonalDetector.Initialize(setup);
			ArrayXf input(17);
			input << 6.101540931246211f, 2.157399825081734f, 0.030358307769100f, 0.001622377620863f, 0.000002650702480f, 0.000082119432769f, 0.000023003709579f, 0.000142849128474f, 0.000095319771848f, 0.000007427525950f, 0.000011033106498f, 0.000002778327015f, 0.000001644832597f, 0.000002744397722f, 0.000003318713830f, 0.000000526210865f, 0.000000215532469f;
			ArrayXf output(17), ref(17);
			tonalDetector.Process(input, output);
			ref << 0.003996608435458f, 0.003674185813906f, 0.001952924840613f, 0.000054817738844f, 0.000014510207955f, 0.000076514186583f, 0.001310586404387f, 0.002839601276986f, 0.002610176308989f, 0.002613706750438f, 0.000543582076583f, 0.001060218869113f, 0.000708952693090f, 0.000726840349593f, 0.002122924273266f, 0.002433458134008f, 0.000517824814604f;
			float error = (ref - output).abs2().sum() / ref.abs2().sum();
			Assert::IsTrue(error < 1e-10f);
		}
	};

	TEST_CLASS(EchoSuppressionCovarianceTest)
	{
		TEST_METHOD(Interface)
		{
			outputLog << "Running EchoSuppressionCovarianceTest->Interface.\n";
			EchoSuppressionCovariance echoSuppression;
			auto c = echoSuppression.GetCoefficients();
			ArrayXXcf input1(c.NBands, c.NChannels);
			ArrayXXcf loopback(c.NBands, c.NChannelsLoopback);
			ArrayXXcf micEstimated(c.NBands, c.NChannels);
			input1.setRandom();
			loopback.setRandom();
			micEstimated.setRandom();
			I::EchoSuppressionCovariance input = { input1, loopback, micEstimated };
			ArrayXXcf output(c.NBands, c.NChannels);
			auto flag = AlgorithmInterfaceTest<EchoSuppressionCovariance>(input, output);
			Assert::IsTrue(flag);
		}

		TEST_METHOD(FreezeLowInput)
		{
			outputLog << "Running EchoSuppressionCovarianceTest->FreezeLowInput.\n";
			EchoSuppressionCovariance echoSuppression;
			auto c = echoSuppression.GetCoefficients();
			c.NChannels = 1;
			echoSuppression.Initialize(c);
			auto p = echoSuppression.GetParameters();
			p.GainReleaseTConstant = p.GainAttackTConstant; // set release gain as fast as attack
			echoSuppression.SetParameters(p);

			ArrayXcf input(c.NBands, c.NChannels);
			input.setConstant(1e-5f);
			ArrayXXcf loopback(c.NBands, c.NChannelsLoopback);
			loopback.setConstant(1e-10f);
			ArrayXXcf micEstimated(c.NBands, c.NChannels);
			micEstimated.setRandom(); // doesnt matter
			ArrayXXcf output(c.NBands, c.NChannels);
			
			int NFrames = 2; // reach Gain=1.f in 2 frames
			for (int i = 0; i < NFrames; i++)
			{
				echoSuppression.Process({ input, loopback, micEstimated }, output);
			}
			float error = (input - output).abs2().sum();
			outputLog << "Error: " << error << std::endl;
			Assert::IsTrue(error < 1e-10f);

			NFrames = 2; //reach minimum gain in 2 frames
			input.setConstant(1e-14f);
			for (int i = 0; i < NFrames; i++)
			{
				echoSuppression.Process({ input, loopback, micEstimated }, output);
			}
			error = (input - output * powf(10.f, -p.GainMinimumdB*.05f)).abs2().sum();
			outputLog << "Error: " << error << std::endl;
			Assert::IsTrue(error < 1e-10f);
		}
	};

	TEST_CLASS(CircBufferTest)
	{
		TEST_METHOD(Interface)
		{
			CircBuffer effect;
			auto c = effect.GetCoefficients();
			int bufferSize = 128;
			ArrayXXf input(bufferSize, c.NChannels);
			input.setRandom();
			ArrayXXf output(bufferSize, c.NChannels);
			auto flag = AlgorithmInterfaceTest<CircBuffer>(input, output);
			Assert::IsTrue(flag);
		}

		TEST_METHOD(ProcessData)
		{
			// process a full delayline 3 times and check that output is equal to input last time
			CircBuffer effect;
			auto c = effect.GetCoefficients();
			effect.Initialize();
			int TimesLength = 3;
			ArrayXXf input(TimesLength * c.DelayLength, c.NChannels);
			input.setRandom();
			ArrayXXf output(c.DelayLength, c.NChannels);
			for (auto i = 0;i < TimesLength;i++)
			{
				effect.Process(input.middleRows(i*c.DelayLength, c.DelayLength), output);
				if (i > 0)
				{
					float error = (output - input.middleRows((i - 1)*c.DelayLength, c.DelayLength)).abs2().sum();
					Assert::IsTrue(error < 1e-10f);
				}
			}
		}

		TEST_METHOD(ReadWriteData)
		{
			// Process some data into delay line and use the different () operators to read the data. Then push one additional value in, and check that the read index has moved correctly
			CircBuffer delay;
			auto c = delay.GetCoefficients();
			c.DelayLength = 4;
			c.NChannels = 2;
			delay.Initialize(c);
			ArrayXXf input(c.DelayLength, c.NChannels);
			input << 0, 4, 1, 5, 2, 6, 3, 7;
			ArrayXXf output(c.DelayLength, c.NChannels);
			delay.Process(input, output);
			//outputLog << delay(3, 1) << ", " << delay(2) << ", " << delay(.7f) << ", " << delay(1.3f, 1) << std::endl;
			Assert::IsTrue(std::abs(delay(3, 1) - 4) < 1e-10f);
			Assert::IsTrue((delay(2) - delay(0.7f) + 1.3f).abs2().sum() < 1e-10f);
			Assert::IsTrue(std::abs(delay(1.3f, 1) - 5.7f) < 1e-10f);

			delay.Push(input.row(0));
			//outputLog << delay(3, 1) << ", " << delay(2) << ", " << delay(.7f) << ", " << delay(1.3f, 1) << std::endl;
			Assert::IsTrue(std::abs(delay(3, 1) - 5) < 1e-10f);
			Assert::IsTrue((delay(2) - delay(0.7f) + .1f).abs2().sum() < 1e-10f);
			Assert::IsTrue(std::abs(delay(1.3f, 1) - 6.7f) < 1e-10f);
		}
	};

	TEST_CLASS(IIR2ndDFTest)
	{
		TEST_METHOD(InterfaceStreaming2ndCascaded)
		{
			IIR2ndCascaded filter;
			auto c = filter.GetCoefficients();
			ArrayXXf sos(c.Nsos, 6);
			sos.setRandom();
			ArrayXXf input(128, c.NChannels), output(128, c.NChannels);
			input.setRandom();
			float sampleRate = 16e3f;
			AlgorithmInterfaceStreamingTest<IIR2ndCascadedStreaming>(input, output, sampleRate, sos);
		}

		TEST_METHOD(CheckCalculationIIR2ndCascaded)
		{
			IIR2ndCascaded filter;
			auto c = filter.GetCoefficients();
			c.NChannels = 1;
			c.Nsos = 3;
			filter.Initialize(c);
			ArrayXXf sos(c.Nsos, 6);
			sos << 1.f, 1.f, -1.f, 1.f, .2f, -.2f, 1.f, -1.f, -1.f, 1.f, .3f, -.4f, 1.f, 2.f, 1.f, 1.f, -.1f, .05f;
			filter.SetPersistentInput(sos);
			ArrayXf input(16), output(16), outputRef(16);
			input.setZero();
			input(0) = 1.f;
			outputRef << 1.000000000000000f, 1.600000000000000f, -2.100000000000000f, -4.234999999999999f, -1.409400000000000f, 0.056560000000000f, 0.118665000000000f, 0.142948000000000f, 0.049654860000000f, 0.056578021000000f, 0.011437689000000f, 0.020982221200000f, -0.000422592759000f, 0.008579439073600f, -0.002496334514100f, 0.004145097177685f;
			filter.Process(input, output);
			float error = (output - outputRef).abs2().sum() / outputRef.abs2().sum();
			outputLog << "error: " << error << std::endl;
			Assert::IsTrue(error < 1e-10f);
		}
		TEST_METHOD(InterfaceStreaming2ndOrder)
		{
			outputLog << "Running IIRFiltersTest->InterfaceStreaming2ndOrder.\n";
			IIR2ndDF2 filter;
			auto c = filter.GetCoefficients();
			auto bufferSize = 128;
			ArrayXXf input(bufferSize, c.NChannels), output(bufferSize, c.NChannels);
			input.setRandom();
			float sampleRate = c.SampleRate;
			auto flag = AlgorithmInterfaceStreamingTest<IIR2ndDF2Streaming>(input, output, sampleRate);
			Assert::IsTrue(flag);

			IIR2ndDF2Transposed filterII;
			c = filterII.GetCoefficients();
			input.resize(bufferSize, c.NChannels);
			output.resize(bufferSize, c.NChannels);
			input.setRandom();
			sampleRate = c.SampleRate;
			flag = AlgorithmInterfaceStreamingTest<IIR2ndDF2TransposedStreaming>(input, output, sampleRate);
			Assert::IsTrue(flag);
		}

		TEST_METHOD(Difference2ndOrder)
		{
			IIR2ndDF2 filterI;
			IIR2ndDF2Transposed filterII;
			auto c = filterI.GetCoefficients();
			filterI.Initialize();
			filterII.Initialize();

			auto pI = filterI.GetParameters();
			auto pII = filterII.GetParameters();
			pI.FilterType = pI.Highpass;
			pII.FilterType = pII.Highpass;
			filterI.SetParameters(pI);
			filterII.SetParameters(pII);
			auto bufferSize = 128;

			ArrayXXf input(bufferSize, c.NChannels), outputI(bufferSize, c.NChannels), outputII(bufferSize, c.NChannels);
			auto NFrames = static_cast<int>(1 * c.SampleRate / bufferSize);
			for (auto i = 0;i < NFrames;i++)
			{
				input.setRandom();
				filterI.Process(input, outputI);
				filterII.Process(input, outputII);
				Assert::IsTrue((input - outputI).abs2().sum() > 1e-10);
				Assert::IsTrue((outputII - outputI).abs2().sum() < 1e-10f);
			}
		}
	};

	TEST_CLASS(IIR2ndNonLinTest)
	{
		TEST_METHOD(InterfaceStreamingNonLinLowpass)
		{
			outputLog << "Running IIRFiltersTest->InterfaceStreamingNonLinLowpass.\n";
			IIR2ndNonLinLowpass filter;
			auto c = filter.GetCoefficients();
			auto bufferSize = 128;
			ArrayXXf input(bufferSize, c.NChannels), output(bufferSize, c.NChannels);
			input.setRandom();
			output.setZero();
			float sampleRate = c.SampleRate;
			auto flag = AlgorithmInterfaceStreamingTest<IIR2ndNonLinLowpassStreaming>(input, output, sampleRate);
			Assert::IsTrue(flag);
		}

		TEST_METHOD(InterfaceStreamingNonLinBandpass)
		{
			outputLog << "Running IIRFiltersTest->InterfaceStreamingNonLinBandpass.\n";
			IIR2ndNonLinBandpass filter;
			auto c = filter.GetCoefficients();
			auto bufferSize = 128;
			ArrayXXf input(bufferSize, c.NChannels), output(bufferSize, c.NChannels);
			input.setRandom();
			output.setZero();
			float sampleRate = c.SampleRate;
			auto flag = AlgorithmInterfaceStreamingTest<IIR2ndNonLinBandpassStreaming>(input, output, sampleRate);
			Assert::IsTrue(flag);
		}

		TEST_METHOD(CalculationNonLinLowpass)
		{
			IIR2ndNonLinLowpass filter;
			auto s = filter.GetSetup();
			int bufferSize = 5;
			float sampleRate = 44100.f;
			s.Coefficients.SampleRate = sampleRate;
			s.Coefficients.NChannels = 1;
			s.Parameters.Resonance = 4.f;
			s.Parameters.Cutoff = 1000.f;
			filter.SetSetup(s);
			filter.Initialize();

			ArrayXf input(bufferSize), output(bufferSize), ref(bufferSize);
			input << 0.1f, 0.2f, 0.3f, 0.4f, 0.5f;
			filter.Process(input, output);
			ref << 0.000241924f, 0.00229658f, 0.00895665f, 0.0237188f, 0.0505466f;
			auto error = (ref - output).abs2().sum() / ref.abs2().sum();
			outputLog << "Error: " << error << std::endl;
			Assert::IsTrue(error < 1e-10f);
		}
	};

	TEST_CLASS(IIR2ndTimeVaryingTest)
	{
		TEST_METHOD(InterfaceStreamingTimeVarying)
		{
			outputLog << "Running IIRFiltersTest->InterfaceTimeStreamingVarying.\n";
			IIR2ndTimeVaryingFilter filter;
			auto c = filter.GetCoefficients();
			auto bufferSize = 128;
			ArrayXXf input(bufferSize, c.NChannels), BandPass(bufferSize, c.NChannels), HighPass(bufferSize, c.NChannels), LowPass(bufferSize, c.NChannels);
			input.setRandom();
			O::IIR2ndTimeVaryingFilter output = { HighPass, LowPass, BandPass };
			auto flag = AlgorithmInterfaceTest<IIR2ndTimeVaryingFilter>(input, output);
			Assert::IsTrue(flag);

			IIR2ndTimeVaryingLowPassFilter filterLowPass;
			auto cLP = filterLowPass.GetCoefficients();
			ArrayXXf outputLP(bufferSize, cLP.NChannels);
			float sampleRate = cLP.SampleRate;
			flag = AlgorithmInterfaceStreamingTest<IIR2ndTimeVaryingLowPassFilterStreaming>(input, outputLP, sampleRate);
			Assert::IsTrue(flag);

			IIR2ndTimeVaryingHighPassFilter filterHighPass;
			auto cHP = filterHighPass.GetCoefficients();
			ArrayXXf outputHP(bufferSize, cHP.NChannels);
			sampleRate = cHP.SampleRate;
			flag = AlgorithmInterfaceStreamingTest<IIR2ndTimeVaryingHighPassFilterStreaming>(input, outputHP, sampleRate);
			Assert::IsTrue(flag);

			IIR2ndTimeVaryingBandPassFilter filterBandPass;
			auto cBP = filterBandPass.GetCoefficients();
			ArrayXXf outputBP(bufferSize, cBP.NChannels);
			sampleRate = cBP.SampleRate;
			flag = AlgorithmInterfaceStreamingTest<IIR2ndTimeVaryingBandPassFilterStreaming>(input, outputBP, sampleRate);
			Assert::IsTrue(flag);

			IIR2ndTimeVaryingNotchFilter filterNotch;
			auto cN = filterNotch.GetCoefficients();
			ArrayXXf outputN(bufferSize, cN.NChannels);
			sampleRate = cN.SampleRate;
			flag = AlgorithmInterfaceStreamingTest<IIR2ndTimeVaryingNotchFilterStreaming>(input, outputN, sampleRate);
			Assert::IsTrue(flag);

			IIR2ndTimeVaryingBellFilter filterBell;
			auto cB = filterBell.GetCoefficients();
			ArrayXXf outputB(bufferSize, cB.NChannels);
			sampleRate = cB.SampleRate;
			flag = AlgorithmInterfaceStreamingTest<IIR2ndTimeVaryingBellFilterStreaming>(input, outputB, sampleRate);
			Assert::IsTrue(flag);

			IIR2ndTimeVaryingHighShelfFilter filterHighShelf;
			auto cHS = filterHighShelf.GetCoefficients();
			ArrayXXf outputHS(bufferSize, cHS.NChannels);
			sampleRate = cHS.SampleRate;
			flag = AlgorithmInterfaceStreamingTest<IIR2ndTimeVaryingHighShelfFilterStreaming>(input, outputHS, sampleRate);
			Assert::IsTrue(flag);

			IIR2ndTimeVaryingLowShelfFilter filterLowShelf;
			auto cLS = filterLowShelf.GetCoefficients();
			ArrayXXf outputLS(bufferSize, cLS.NChannels);
			sampleRate = cLS.SampleRate;
			flag = AlgorithmInterfaceStreamingTest<IIR2ndTimeVaryingLowShelfFilterStreaming>(input, outputLS, sampleRate);
			Assert::IsTrue(flag);
		}
	};

	TEST_CLASS(FFTTest)
	{
		TEST_METHOD(Interface)
		{
			FFTReal FFT;
			auto c = FFT.GetCoefficients();
			auto p = FFT.GetParameters();
			auto nChannels = 2;
			ArrayXXf input(c.FFTSize, nChannels);
			input.setRandom();
			ArrayXXcf output(c.FFTSize / 2 + 1, nChannels);
			auto flag = AlgorithmInterfaceTest<FFTReal>(input, output);
			Assert::IsTrue(flag);
		}

		TEST_METHOD(InterfaceInverse)
		{
			FFTRealInverse FFTInverse;
			auto c = FFTInverse.GetCoefficients();
			auto p = FFTInverse.GetParameters();
			auto nChannels = 2;
			ArrayXXcf input(c.FFTSize/2+1, nChannels);
			input.setRandom();
			ArrayXXf output(c.FFTSize, nChannels);
			auto flag = AlgorithmInterfaceTest<FFTRealInverse>(input, output);
			Assert::IsTrue(flag);
		}

		TEST_METHOD(CompareInverse)
		{
			FFTReal FFT;
			FFTRealInverse FFTInverse;
			auto s = FFT.GetSetup();
			auto sInv = FFTInverse.GetSetup();
			sInv.Coefficients.FFTSize = s.Coefficients.FFTSize;
			FFT.Initialize(s);
			FFTInverse.Initialize(sInv);

			ArrayXcf input(s.Coefficients.FFTSize / 2 + 1);
			input.setRandom();
			input.imag()(0) = 0.f;
			input.imag()(s.Coefficients.FFTSize/2) = 0.f;
			ArrayXf output(s.Coefficients.FFTSize), output2(s.Coefficients.FFTSize);
			FFT.Inverse(input, output);
			FFTInverse.Process(input, output2);
			float error = (output - output2).abs2().sum();
			outputLog << "Error: " << error << std::endl;
			Assert::IsTrue(error == 0.f);
		}

		TEST_METHOD(ForwardBackward)
		{
			FFTReal FFT;
			FFTRealInverse FFTInverse;
			int fftSize = 32;
			auto s = FFT.GetSetup();
			auto sInv = FFTInverse.GetSetup();
			s.Coefficients.FFTSize = fftSize;
			sInv.Coefficients.FFTSize = fftSize;
			FFT.Initialize(s);
			FFTInverse.Initialize(sInv);

			ArrayXf input(s.Coefficients.FFTSize), res(s.Coefficients.FFTSize);
			input.setRandom();
			ArrayXcf output(s.Coefficients.FFTSize/2+1);
			FFT.Process(input, output);
			FFTInverse.Process(output, res);
			float error = (input-res).abs2().sum();
			outputLog << "Error: " << error << std::endl;
			Assert::IsTrue(error < 1e-10f);
		}
	};

	TEST_CLASS(FFTConvolutionTest)
	{
		TEST_METHOD(InterfaceStreaming)
		{
			FFTConvolution conv;
			auto c = conv.GetCoefficients();
			ArrayXXf input(c.BufferSize, c.NChannels);
			ArrayXf filter(c.FilterSize);
			input.setRandom();
			filter.setRandom();
			ArrayXXf output(c.BufferSize, c.NChannels * c.NFiltersPerChannel);
			float sampleRate = 48000.f;
			auto flag = AlgorithmInterfaceStreamingTest<FFTConvolutionStreaming>(input, output, sampleRate, filter);
			Assert::IsTrue(flag);
		}
	};

	TEST_CLASS(FilterbankTest)
	{
		TEST_METHOD(InterfaceAnalysis)
		{
			outputLog << "Running FilterbankTest->InterfaceAnalysis.\n";
			FilterbankAnalysis filterbank;
			auto c = filterbank.GetCoefficients();
			ArrayXXf input(c.BufferSize, c.NChannels);
			input.setRandom();
			ArrayXXcf output(c.FFTSize / 2 + 1, c.NChannels);
			Assert::IsTrue(AlgorithmInterfaceTest<FilterbankAnalysis>(input, output));
		}

		TEST_METHOD(InterfaceSynthesis)
		{
			outputLog << "Running FilterbankTest->InterfaceSynthesisPFFFT.\n";
			FilterbankSynthesis filterbank;
			auto c = filterbank.GetCoefficients();
			ArrayXXcf input(c.FrameSize / 2 + 1, c.NChannels);
			input.setRandom();
			ArrayXXf output(c.BufferSize, c.NChannels);
			Assert::IsTrue(AlgorithmInterfaceTest<FilterbankSynthesis>(input, output));
		}

		TEST_METHOD(PassThrough)
		{
			outputLog << "Running FilterbankWindowTest->PassThrough.\n";
			FilterbankAnalysis filterbank;
			FilterbankSynthesis filterbankInverse;
			auto cA = filterbank.GetCoefficients();
			auto cS = filterbankInverse.GetCoefficients();
			cS.NChannels = cA.NChannels;
			filterbank.Initialize(cA);
			filterbankInverse.Initialize(cS);
			const int NFrames = 30;
			ArrayXXf input(NFrames*cA.BufferSize, cA.NChannels);
			input.setRandom();
			//input.row(0) = Eigen::ArrayXf::LinSpaced(NFrames*cA.BufferSize, 0, NFrames*cA.BufferSize - 1);
			ArrayXXf output(NFrames*cA.BufferSize, cS.NChannels);
			ArrayXXcf xFreq(cA.FrameSize / 2 + 1, cA.NChannels);
			for (auto i = 0; i < NFrames; i++)
			{
				filterbank.Process(input.middleRows(i*cA.BufferSize, cA.BufferSize), xFreq);
				filterbankInverse.Process(xFreq, output.middleRows(i*cS.BufferSize, cS.BufferSize));
			}
			auto error = (output.bottomRows(NFrames*cA.BufferSize - cA.FrameSize + cA.BufferSize) - input.topRows(NFrames*cA.BufferSize - cA.FrameSize + cA.BufferSize)).abs2().mean();
			outputLog << "Error: " << error << std::endl;
			Assert::IsTrue(error < 1e-10f);
			outputLog << "Test succesful." << std::endl;
		}
	};

	TEST_CLASS(FIRMinPhaseTest)
	{
		TEST_METHOD(Interface)
		{
			outputLog << "Running FIRMinPhaseTest->Interface.\n";
			FIRMinPhase filter;
			auto c = filter.GetCoefficients();
			ArrayXf input(c.NChannels);
			input.setRandom();
			ArrayXf output(c.NChannels);
			auto flag = AlgorithmInterfaceTest<FIRMinPhase>(input, output);
			Assert::IsTrue(flag);
		}
	};

	TEST_CLASS(EchoCancellerMomentumTest)
	{
		class EchoCancellerMomentumDebug : public EchoCancellerMomentum
		{
		public:
			using EchoCancellerMomentum::GetPrivateData;
		};
		TEST_METHOD(Interface)
		{
			outputLog << "Running EchoCancellerMomentumTest->Interface.\n";
			EchoCancellerMomentum echoCanceller;
			auto c = echoCanceller.GetCoefficients();
			ArrayXXcf input1(c.NBands, c.NChannels);
			ArrayXXcf input2(c.NBands, c.NChannelsLoopback);
			input1.setRandom();
			input2.setRandom();
			I::EchoCancellerMomentum input = { input1, input2 };
			ArrayXXcf output(c.NBands, c.NChannels);
			Assert::IsTrue(AlgorithmInterfaceTest<EchoCancellerMomentum>(input, output));
		}

		TEST_METHOD(FilterSubbands)
		{
			outputLog << "Running EchoCancellerMomentumTest->Filter.\n";
			EchoCancellerMomentumDebug echoCanceller;
			auto c = echoCanceller.GetCoefficients();
			c.NChannelsLoopback = 1;
			c.NChannels = 1;
			c.NBands = 10;
			echoCanceller.Initialize(c);
			auto P = echoCanceller.GetParameters();
			P.FreezeUpdate = true;
			echoCanceller.SetParameters(P);
			echoCanceller.EchoSuppression.Disable();
			auto Data = echoCanceller.GetPrivateData();
			Data->Filters[0][0].setZero();
			Data->Filters[0][0].col(1) = 1.f;
			Data->Filters[0][0].col(2) = 2.f;
			int NFrames = Data->BufferLength * 10;
			ArrayXXcf Input(c.NBands * 2, NFrames), Loopback(c.NBands * 2, NFrames), Output(c.NBands * 2, NFrames);

			Map<ArrayXXcf> InputSubbands = Map<ArrayXXcf>(nullptr, 0, 0);
			Map<ArrayXXcf> LoopbackSubbands = Map<ArrayXXcf>(nullptr, 0, 0);
			Map<ArrayXXcf> OutputSubbands = Map<ArrayXXcf>(nullptr, 0, 0);

			new (&InputSubbands) Eigen::Map<Eigen::ArrayXXcf>(&Input(0, 0), c.NBands, NFrames);
			new (&LoopbackSubbands) Eigen::Map<Eigen::ArrayXXcf>(&Loopback(0, 0), c.NBands, NFrames);
			new (&OutputSubbands) Eigen::Map<Eigen::ArrayXXcf>(&Output(0, 0), c.NBands, NFrames);

			LoopbackSubbands.setRandom();
			InputSubbands.setZero();
			OutputSubbands.setZero();
			for (auto t = 0; t < NFrames; t++)
			{
				echoCanceller.Process({ InputSubbands.col(t), LoopbackSubbands.col(t) }, OutputSubbands.col(t));
			}
			ArrayXXcf Res = -LoopbackSubbands.leftCols(LoopbackSubbands.cols() - 1);
			Res.rightCols(Res.cols() - 1) -= 2 * LoopbackSubbands.leftCols(LoopbackSubbands.cols() - 2);
			float error = (OutputSubbands.rightCols(OutputSubbands.cols() - 1) - Res).abs2().sum();
			outputLog << "Error: " << error << std::endl;
			Assert::IsTrue(error == 0);
			outputLog << "Test succesful." << std::endl;
		}

		TEST_METHOD(LessThanTwoLoopbackChannels)
		{
			outputLog << "Running EchoCancellerMomentumTest->LessThanTwoLoopbackChannels.\n";
			EchoCancellerMomentum echoCanceller;
			auto c = echoCanceller.GetCoefficients();
			ArrayXXcf input(c.NBands, c.NChannels), output(c.NBands, c.NChannels);
			input.setRandom();

			c.NChannelsLoopback = 1;
			echoCanceller.Initialize(c);
			echoCanceller.EchoSuppression.Disable();
			ArrayXXcf Loopback(c.NBands, c.NChannelsLoopback);
			Loopback.setRandom();
			echoCanceller.Process({ input, Loopback }, output);

			c.NChannelsLoopback = 0;
			echoCanceller.Initialize(c);
			echoCanceller.EchoSuppression.Disable();
			Loopback.resize(c.NBands, c.NChannelsLoopback);
			Loopback.setRandom();
			echoCanceller.Process({ input, Loopback }, output);
			float error = (output - input).abs2().sum();
			Assert::IsTrue(error == 0.f);
			outputLog << "Test succesful." << std::endl;
		}

		TEST_METHOD(Process)
		{
			using namespace std::literals;

			outputLog << "Running EchoCancellerMomentumTest->Process.\n";
			srand(1);
			EchoCancellerMomentumDebug echoCanceller;
			auto c = echoCanceller.GetCoefficients();
			c.NChannels = 1;
			c.NChannelsLoopback	= 2;
			echoCanceller.Initialize(c);
			auto P = echoCanceller.GetParameters();
			P.FreezeUpdate = true;
			echoCanceller.SetParameters(P);
			echoCanceller.EchoSuppression.Disable();
			auto Data = echoCanceller.GetPrivateData();
			auto BufferLength = static_cast<int>(Data->Filters[0][0].cols());
			for (auto &filter : Data->Filters[0]) { filter = ArrayXXcf::Random(c.NBands, BufferLength) + 1i*ArrayXXcf::Random(c.NBands, BufferLength); }
			ArrayXXcf loopback(c.NBands, c.NChannelsLoopback);

			int NFrames = 20;
			std::vector<ArrayXXcf> loopbackData(2);
			for (auto &loopback : loopbackData)
			{
				loopback = ArrayXXcf::Random(c.NBands, NFrames) + 1i*ArrayXXcf::Random(c.NBands, NFrames);
			}
			ArrayXXcf input = ArrayXXcf::Random(c.NBands, NFrames) + 1i*ArrayXXcf::Random(c.NBands, NFrames);
			ArrayXXcf output = ArrayXXcf::Zero(c.NBands, NFrames);
			for (int i = 0; i < NFrames; i++)
			{
				loopback.col(0) = loopbackData[0].col(i);
				loopback.col(1) = loopbackData[1].col(i);
				echoCanceller.Process({ input.col(i), loopback }, output.col(i));

				ArrayXcf temp = input.col(i) - output.col(i);
				ArrayXXcf temp2 = loopbackData[0].middleCols((std::max)(0, -BufferLength + i + 1), (std::min)(i + 1, BufferLength))*Data->Filters[0][0].leftCols((std::min)(i + 1, BufferLength)).rowwise().reverse();
				temp2 += loopbackData[1].middleCols((std::max)(0, -BufferLength + i + 1), (std::min)(i + 1, BufferLength))*Data->Filters[0][1].leftCols((std::min)(i + 1, BufferLength)).rowwise().reverse();
				float error = (temp - temp2.rowwise().sum()).abs2().sum();
				Assert::IsTrue(error < 1e-8f);
			}
			outputLog << "Test succesful." << std::endl;
		}

#if NDEBUG 
		TEST_METHOD(Update) // this test calculates differently in release and debug mode
		{
			outputLog << "Running EchoCancellerMomentumTest->Update.\n";
			EchoCancellerMomentumDebug echoCanceller;
			auto c = echoCanceller.GetCoefficients();
			c.NChannels = 1;
			c.NBands = 10; // reduce number of bands just to make calculations more efficient
			c.FilterbankRate = 100.f;
			echoCanceller.Initialize(c);
			echoCanceller.EchoSuppression.Disable();
			auto Data = echoCanceller.GetPrivateData();

			ArrayXcf input(c.NBands);
			ArrayXcf output(c.NBands);
			ArrayXXcf loopback(c.NBands, 2);

			int NFrames = 5;
			for (int i = 0; i < NFrames; i++)
			{
				input.setConstant(1.f - i * .25f);
				loopback.col(0).setConstant(i*.25f);
				loopback.col(1).setConstant(1.f - i * .25f);
				echoCanceller.Process({ input, loopback }, output);
			}
			ArrayXcf true0(4), true1(4);
			true0 << 0.001110948529850f, 0.000450180134629f, 0.000112396285537f, -0.000001563287774f;
			true1 << 0.003180645515587f, 0.002190628969350f, 0.001237722373967f, 0.000457482862416f;
			float error = (Data->Filters[0][0].row(0).transpose() - true0).abs2().sum();
			error += (Data->Filters[0][1].row(0).transpose() - true1).abs2().sum();
			Assert::IsTrue(error < 1e-10);
			outputLog << "Test succesful." << std::endl;
		}
#endif // #if NDEBUG
	};

	TEST_CLASS(Upsampling2XCubicTest)
	{
		TEST_METHOD(Interface2x)
		{
			outputLog << "Running Upsampling2XCubicTest->Interface.\n";
			Upsampling2XCubic upsampling;
			ArrayXXf input(256, 2), output(512, 2);
			input.setRandom();
			auto flag = AlgorithmInterfaceTest<Upsampling2XCubic>(input, output);
			Assert::IsTrue(flag);
		}

		TEST_METHOD(InterfacePower2)
		{
			outputLog << "Running Upsampling2XCubicTest->Interface.\n";
			UpsamplingPower2Cubic upsampling;
			auto p = upsampling.GetParameters();
			ArrayXXf input(256, 2), output(256 * (2 << (p.UpsamplingFactorPow2-1)), 2);
			input.setRandom();
			auto flag = AlgorithmInterfaceTest<UpsamplingPower2Cubic>(input, output);
			Assert::IsTrue(flag);
		}
	};

	TEST_CLASS(VirtualizationHeadphonesTest)
	{
		TEST_METHOD(InterfaceStreaming)
		{
			outputLog << "Running VirtualizationHeadphonesTest->Interface.\n";
			VirtualizationHeadphones virtualHeadphones;
			auto c = virtualHeadphones.GetCoefficients();
			ArrayXXf input(c.BufferSize, c.NChannels), output(c.BufferSize, c.NChannels);
			input.setRandom();
			ArrayXf sub(c.BufferSize);
			sub.setRandom();
			float sampleRate;
			if (c.SampleRate == c.Hz44100) { sampleRate = 44100.f; }
			else if (c.SampleRate == c.Hz48000) { sampleRate = 48000.f; }
			auto flag = AlgorithmInterfaceStreamingTest<VirtualizationHeadphonesStreaming>(input, output, sampleRate);
			Assert::IsTrue(flag);
		}

		TEST_METHOD(CheckImpulseResponseAtDifferentSampleRates)
		{
			outputLog << "Running VirtualizationHeadphonesTest->CheckImpulseResponseAtDifferentSampleRates.\n";

			// process at 48kHz
			VirtualizationHeadphones virtualHeadphonesx1;
			auto c = virtualHeadphonesx1.GetCoefficients();
			c.SampleRate = c.Hz48000;
			virtualHeadphonesx1.Initialize(c);
			ArrayXXf inputx1(c.BufferSize, c.NChannels), output(c.BufferSize, 2);
			inputx1.setZero();
			inputx1(0, 0) = 1.f;
			ArrayXf sub(c.BufferSize);
			sub.setRandom();
			virtualHeadphonesx1.Process(inputx1, output);

			// process at 96kHz
			VirtualizationHeadphones virtualHeadphonesx2;
			auto cx2 = virtualHeadphonesx2.GetCoefficients();
			cx2.SampleRate = cx2.Hz96000;
			cx2.BufferSize = 2 * cx2.BufferSize;
			virtualHeadphonesx2.Initialize(cx2);
			ArrayXXf inputx2(cx2.BufferSize, cx2.NChannels), outputx2(cx2.BufferSize, 2);
			inputx2.setZero();
			inputx2(0, 0) = 1.f;
			ArrayXf subx2(cx2.BufferSize);
			subx2.setRandom();
			virtualHeadphonesx2.Process(inputx2, outputx2);
			
			float error = 0;
			for (auto i = 0; i < 11; i++) // only check before the highpass filtered impulse is added to the 96kHz impulse response
			{
				error += std::abs(outputx2(2 * i,0)*2 - output(i,0));
			}
			outputLog << "Error: " << error << std::endl;
			Assert::IsTrue(error < 1e-7f);

			// process at 192kHz
			VirtualizationHeadphones virtualHeadphonesx4;
			auto cx4 = virtualHeadphonesx4.GetCoefficients();
			cx4.SampleRate = cx4.Hz192000;
			virtualHeadphonesx4.Initialize(cx4);
			ArrayXXf inputx4(cx4.BufferSize, cx4.NChannels), outputx4(cx4.BufferSize, 2);
			inputx4.setZero();
			inputx4(0, 0) = 1.f;
			ArrayXf subx4(cx4.BufferSize);
			subx4.setRandom();
			virtualHeadphonesx4.Process(inputx4, outputx4);

			error = 0;
			for (auto i = 0; i < 10; i++) // only check before the highpass filtered impulse is added to the 192kHz impulse response
			{
				error += std::abs(outputx4(4 * i,0) * 4 - output(i,0));
			}
			outputLog << "Error: " << error << std::endl;
			Assert::IsTrue(error < 1e-7f);

			// process at 44100kHz
			c.SampleRate = c.Hz44100;
			virtualHeadphonesx1.Initialize(c);
			virtualHeadphonesx1.Process(inputx1, output);

			// process at 88200kHz
			cx2.SampleRate = cx2.Hz88200;
			virtualHeadphonesx2.Initialize(cx2);
			virtualHeadphonesx2.Process(inputx2, outputx2);

			error = 0;
			for (auto i = 0; i < 11; i++) // only check before the highpass filtered impulse is added to the 88.2kHz impulse response
			{
				error += std::abs(outputx2(2 * i, 0) * 2 - output(i, 0));
			}
			outputLog << "Error: " << error << std::endl;
			Assert::IsTrue(error < 1e-7f);

			// process at 176400kHz
			cx4.SampleRate = cx4.Hz176400;
			virtualHeadphonesx4.Initialize(cx4);
			virtualHeadphonesx4.Process(inputx4, outputx4);

			error = 0;
			for (auto i = 0; i < 10; i++) // only check before the highpass filtered impulse is added to the 176.4kHz impulse response
			{
				error += std::abs(outputx4(4 * i, 0) * 4 - output(i, 0));
			}
			outputLog << "Error: " << error << std::endl;
			Assert::IsTrue(error < 1e-7f);
		}

		TEST_METHOD(CheckSub)
		{
			outputLog << "Running VirtualizationHeadphonesTest->CheckSub.\n";

			VirtualizationHeadphones virtualHeadphones;
			auto c = virtualHeadphones.GetCoefficients();
			c.EnabledSub = true;
			virtualHeadphones.Initialize(c);
			int delay = virtualHeadphones.GetLatencySamples();

			ArrayXXf input(c.BufferSize, c.NChannels), output(c.BufferSize,2);
			input.setZero();
			virtualHeadphones.Process(input, output);
			float error = (output.bottomRows(c.BufferSize - delay) - input.rightCols<1>().head(c.BufferSize - delay).replicate(1,2)).abs2().sum();
			outputLog << "Error: " << error << std::endl;
			Assert::IsTrue(error < 1e-10f);

			virtualHeadphones.Process(input, output);
			error = (output.topRows(delay) - input.rightCols<1>().tail(delay).replicate(1, 2)).abs2().sum();
			outputLog << "Error: " << error << "\n";
			error += (output.bottomRows(c.BufferSize - delay) - input.rightCols<1>().head(c.BufferSize - delay).replicate(1, 2)).abs2().sum();
			outputLog << "Error: " << error << std::endl;
			Assert::IsTrue(error < 1e-10f);
		}
	};

	TEST_CLASS(VoiceActivityDetectionTest)
	{

		TEST_METHOD(Interface)
		{
			outputLog << "Running VoiceActivityDetectionTest->Interface.\n";
			VoiceActivationDetection VAD;
			auto c = VAD.GetCoefficients();
			ArrayXXcf input(c.NBands, c.NChannels);
			auto output = true;
			input.setRandom();
			auto flag = AlgorithmInterfaceTest<VoiceActivationDetection>(input, output);
			Assert::IsTrue(flag);
		}

		TEST_METHOD(ActivityDetect)
		{
			outputLog << "Running VoiceActivityDetectionTest->ActivityDetect.\n";
			VoiceActivationDetection VAD;
			auto c = VAD.GetCoefficients();
			c.NBands = 256;
			VAD.Initialize(c);

			// run for .8 seconds with random noise
			int Nend = int(.8f * c.FilterbankRate);
			bool Activity = false;
			for (auto n = 0; n < Nend; n++)
			{
				ArrayXXcf Noise = 1e-2f * ArrayXXcf::Random(c.NBands, c.NChannels);
				VAD.Process(Noise, Activity);
			}

			// Add activity at certain times and test activity detection
			for (auto n = 0; n < 51; n++)
			{
				ArrayXXcf Noise = 1e-2f * ArrayXXcf::Random(c.NBands, c.NChannels);
				switch (n)
				{
				case 0: Noise.col(0) = ArrayXcf::Ones(c.NBands); break;
				case 10:Noise.col(1) = ArrayXcf::Ones(c.NBands); break;
				case 20:Noise.col(0).segment(20, 45) = ArrayXcf::Constant(45, 2.f); break;
				case 30:Noise.col(1).segment(70, 45) = ArrayXcf::Constant(45, 2.f); break;
				case 40:Noise.col(0).segment(133, 45) = ArrayXcf::Constant(45, 2.f); break;
				case 50:Noise.col(1).segment(201, 45) = ArrayXcf::Constant(45, 2.f); break;
				default: break;
				}
				VAD.Process(Noise, Activity);

				if (n % 10 == 0)
				{
					Assert::IsTrue(Activity);
				}
				else {
					Assert::IsFalse(Activity);
				}
			}
			outputLog << "Test succesful." << std::endl;
		}
	};
}

namespace UtilitiesUnitTests
{
	TEST_CLASS(ApproximationMathTest)
	{
	public:

		TEST_METHOD(PowerValueTest)
		{
			outputLog << "Running FunctionsTest->ValueTest: Test accuracy of coarse power calculations for 1000 random numbers.\n";
			for (auto i = 0; i < 1000; i++)
			{
				// Uniform random number between -9 and 9 (10^9 ~ 2^31)
				// The approximations can only correctly calculate values between +-2^31, due to integer precision of exponent
				float x = 18.f * static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX) - 9.f;

				float y2Int = pow(2.f, static_cast<int>(x));
				float y2IntTest = Pow2(static_cast<int>(x));

				float ye = powf(2.718281828459046f, x);
				float yeTest = PowECoarse(x);

				float y10 = powf(10.f, x);
				float y10Test = Pow10Coarse(x);

				float y2 = powf(2.f, x);
				float y2Test = Pow2Coarse(x);

				float yLin = powf(10.f, x / 20.f);
				float yLinTest = Db2LinCoarse(x);

				//char buffer[50];
				//sprintf(buffer, "x: %g, yTrue: %g, yTest: %g", x, y2Int, y2IntTest);
				//Logger::WriteMessage(buffer);
				/*sprintf(buffer, "x: %g, yTrue: %g, yTest: %g", x, ye, yeTest);
				Logger::WriteMessage(buffer);
				sprintf(buffer, "x: %g, yTrue: %g, yTest: %g", x, y10, y10Test);
				Logger::WriteMessage(buffer);
				sprintf(buffer, "x: %g, yTrue: %g, yTest: %g", x, y2, y2Test);
				Logger::WriteMessage(buffer);
				sprintf(buffer, "x: %g, yTrue: %g, yTest: %g", x, yLin, yLinTest);
				Logger::WriteMessage(buffer);*/

				float error = std::abs(y2Int - y2IntTest) / y2Int;
				Assert::IsTrue(error == 0.f);
				error = std::abs(ye - yeTest) / ye;
				Assert::IsTrue(error < 0.0615f);
				error = std::abs(y10 - y10Test) / y10;
				Assert::IsTrue(error < 0.0615f);
				error = std::abs(y2 - y2Test) / y2;
				Assert::IsTrue(error < 0.0615f);
				error = std::abs(yLin - yLinTest) / yLin;
				Assert::IsTrue(error < 0.0615f);
			}
			outputLog << "Test succesful." << std::endl;
		}

		TEST_METHOD(PowerSpeedTest)
		{
			outputLog << "Running FunctionsTest->PowerSpeedTest: Test speed of coarse power calculations.\n";
			int NumberOfIterations = 1001;

			float y = 0.f;
			auto start = std::chrono::high_resolution_clock::now();
			for (auto i = 0; i < NumberOfIterations; i++)
			{
				float x = 2.f * static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX) - 1.f;
				y += powf(10.f, x);
			}
			auto end = std::chrono::high_resolution_clock::now();
			auto diff = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

			outputLog << "pow(10.f,x). Time: " << diff << ", Value: " << y << ".\n";

			y = 0.f;
			start = std::chrono::high_resolution_clock::now();
			for (auto i = 0; i < NumberOfIterations; i++)
			{
				float x = 2.f * static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX) - 1.f;
				y += Pow10Coarse(x);
			}
			end = std::chrono::high_resolution_clock::now();
			diff = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
			outputLog << "Pow10Coarse(10.f,x). Time: " << diff << ", Value: " << y << ".\n";

			outputLog << "Test succesful." << std::endl;
		}

		TEST_METHOD(AddVersusMultTest)
		{
			outputLog << "Running FunctionsTest->AddVersusMult: Test speed of additions versus multiplications.\n";
			int NumberOfIterations = 1001;

			float y = 1.f;
			auto start = std::chrono::high_resolution_clock::now();
			for (auto i = 0; i < NumberOfIterations; i++)
			{
				float x = 4.f * static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX) - 2.f;
				y += x;
			}
			auto end = std::chrono::high_resolution_clock::now();
			auto diff = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

			outputLog << "Summation. Time: " << diff << ", Value:" << y << ".\n";

			y = 1.f;
			start = std::chrono::high_resolution_clock::now();
			for (auto i = 0; i < NumberOfIterations; i++)
			{
				float x = 4.f * static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX) - 2.f;
				y *= x;
			}
			end = std::chrono::high_resolution_clock::now();
			diff = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

			outputLog << "Multiplication. Time: " << diff << ", Value: " << y << ".\n";
			outputLog << "Test succesful." << std::endl;
		}
	};

	TEST_CLASS(FormatHandlingTest)
	{
	public:
		TEST_METHOD(ReadIn)
		{
			outputLog << "Running FormatHandlingTest->ReadIn.\n";
			FormatHandling<float> IOHandler;
			const int Channels = 5;
			const int BufferSize = 100;
			ArrayXXf Buffer;
			Buffer.resize(BufferSize, Channels);
			float Input[Channels * BufferSize];
			for (auto i = 0; i < Channels * BufferSize; i++) { Input[i] = static_cast<float>(i); }

			auto start = std::chrono::steady_clock::now();
			for (int i = 0; i < 100; i++)
			{
				IOHandler.ReadIn(Input, Buffer);
			}

			auto end = std::chrono::steady_clock::now();
			outputLog << "Elapsed Time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << " nanoseconds\n";

			ArrayXXf BufferTest;
			BufferTest.resize(BufferSize, Channels);
			for (auto i = 0; i < Channels; i++)
			{
				BufferTest.col(i) = ArrayXf::LinSpaced(BufferSize, (float)i, (float)(i + (BufferSize - 1) * Channels));
			}
			float error = (Buffer - BufferTest).abs2().sum();
			Assert::IsTrue(error < 1e-10f);

			ArrayXXf BufferRef;
			BufferRef.resize(Channels, BufferSize);
			BufferRef.setRandom();
			start = std::chrono::steady_clock::now();
			for (int i = 0; i < 100; i++)
			{
				BufferRef = Buffer.transpose();
			}
			end = std::chrono::steady_clock::now();
			outputLog << "Elapsed Reference Time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << " nanoseconds\n";

			outputLog << "Test succesful." << std::endl;
		}

		TEST_METHOD(ReadWithStride)
		{
			outputLog << "Running FormatHandlingTest->ReadWithStride.\n";
			FormatHandling<float> IOHandler;
			const int Channels = 5; // stride is 5
			const int ChannelsToRead = 3; // but read channels is only 3
			const int BufferSize = 100;
			ArrayXXf Buffer;
			Buffer.resize(BufferSize, ChannelsToRead);
			float Input[Channels * BufferSize];
			for (auto i = 0; i < Channels * BufferSize; i++) { Input[i] = static_cast<float>(i); }

			auto start = std::chrono::steady_clock::now();
			for (int i = 0; i < 100; i++)
			{
				IOHandler.ReadIn(Channels, Input, Buffer);
			}

			auto end = std::chrono::steady_clock::now();
			outputLog << "Elapsed Time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << " nanoseconds\n";

			ArrayXXf BufferTest;
			BufferTest.resize(BufferSize, ChannelsToRead);
			for (auto i = 0; i < ChannelsToRead; i++)
			{
				BufferTest.col(i) = ArrayXf::LinSpaced(BufferSize, (float)i, (float)(i + (BufferSize - 1) * Channels));
			}
			float error = (Buffer - BufferTest).abs2().sum();
			Assert::IsTrue(error < 1e-10f);

			ArrayXXf BufferRef;
			BufferRef.resize(Channels, BufferSize);
			BufferRef.setRandom();
			start = std::chrono::steady_clock::now();
			for (int i = 0; i < 100; i++)
			{
				BufferRef = Buffer.transpose();
			}
			end = std::chrono::steady_clock::now();
			outputLog << "Elapsed Reference Time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << " nanoseconds\n";

			outputLog << "Test succesful." << std::endl;
		}

		TEST_METHOD(ReadWithLoopback)
		{
			outputLog << "Running FormatHandlingTest->ReadWithLoopback.\n";
			FormatHandling<float> IOHandler;
			const int Channels = 3;
			const int LoopbackChannels = 2;

			const int BufferSize = 100;
			ArrayXXf Buffer;
			Buffer.resize(BufferSize, Channels + LoopbackChannels);
			float Input[Channels * BufferSize];
			float Loopback[LoopbackChannels * BufferSize];
			for (auto i = 0; i < Channels * BufferSize; i++) { Input[i] = static_cast<float>(i); }
			for (auto i = 0; i < LoopbackChannels * BufferSize; i++) { Loopback[i] = static_cast<float>(-i); }

			IOHandler.ReadIn(Input, Buffer.leftCols(Channels));
			IOHandler.ReadIn(Loopback, Buffer.rightCols(LoopbackChannels));

			ArrayXXf BufferTest;
			BufferTest.resize(BufferSize, Channels + LoopbackChannels);
			for (auto i = 0; i < Channels; i++)
			{
				BufferTest.col(i) = ArrayXf::LinSpaced(BufferSize, (float)i, (float)(i + (BufferSize - 1) * Channels));
			}
			for (auto i = 0; i < LoopbackChannels; i++)
			{
				BufferTest.col(i + Channels) = ArrayXf::LinSpaced(BufferSize, (float)(-i), (float)(-i - (BufferSize - 1) * LoopbackChannels));
			}
			float error = (Buffer - BufferTest).abs2().sum();
			Assert::IsTrue(error < 1e-10f);
			outputLog << "Test succesful." << std::endl;
		}

		TEST_METHOD(ReadWithIndex)
		{
			outputLog << "Running FormatHandlingTest->ReadWithIndex.\n";
			FormatHandling<float> IOHandler;
			const int Channels = 5;
			const int BufferSize = 100;
			const int Index = 10;
			ArrayXXf Buffer;
			Buffer.resize(BufferSize, Channels);
			float Input[Channels * (BufferSize + Index)];
			for (auto i = 0; i < Channels * (BufferSize + Index); i++) { Input[i] = static_cast<float>(i); }

			auto start = std::chrono::steady_clock::now();
			IOHandler.ReadIn(Input, Index, Buffer);
			auto end = std::chrono::steady_clock::now();
			outputLog << "Elapsed Time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << " microseconds\n";

			ArrayXXf BufferTest;
			BufferTest.resize(BufferSize, Channels);
			for (auto i = 0; i < Channels; i++)
			{
				BufferTest.col(i) = ArrayXf::LinSpaced(BufferSize, (float)i + Index * Channels, (float)(i + (BufferSize + Index - 1) * Channels));
			}
			float error = (Buffer - BufferTest).abs2().sum();
			Assert::IsTrue(error < 1e-10f);
			outputLog << "Test succesful." << std::endl;
		}

		TEST_METHOD(WriteOut)
		{
			outputLog << "Running FormatHandlingTest->WriteOut.\n";
			FormatHandling<float> IOHandler;
			const int Channels = 5;
			const int BufferSize = 100;
			ArrayXXf Buffer;
			Buffer.resize(BufferSize, Channels);
			for (auto i = 0; i < Channels; i++)
			{
				Buffer.col(i) = ArrayXf::LinSpaced(BufferSize, (float)i, (float)(i + (BufferSize - 1) * Channels));
			}
			float Output[BufferSize];
			float OutputTest[BufferSize];
			for (auto i = 0; i < BufferSize; i++) { OutputTest[i] = static_cast<float>(i*Channels); }

			IOHandler.WriteOut(Buffer.col(0), Output);

			float error = 0;
			for (auto i = 0; i < BufferSize; i++)
			{
				error += fabs(Output[i] - OutputTest[i]);
			}
			Assert::IsTrue(error < 1e-10f);

			outputLog << "Test succesful." << std::endl;
		}

		TEST_METHOD(ReadInInt)
		{
			outputLog << "Running FormatHandlingTest->ReadInInt.\n";
			FormatHandling<int32_t> IOHandler32;
			FormatHandling<int16_t> IOHandler16;
			const int Channels = 4; // Channels must be factor of 2 for casting to int16 doesnt cause quantization
			const int BufferSize = 128; // BufferSize must be factor of 2 for casting to int16 doesnt cause quantization
			ArrayXXf Buffer16, Buffer32;
			Buffer16.resize(BufferSize, Channels);
			Buffer32.resize(BufferSize, Channels);

			int16_t Input16[Channels * BufferSize];
			int32_t Input32[Channels * BufferSize];
			for (auto i = 0; i < Channels * BufferSize; i++)
			{
				Input16[i] = static_cast<int16_t>(i* 32768.f / (Channels * BufferSize));
				Input32[i] = static_cast<int32_t>(Input16[i] * 65536);
			}

			IOHandler16.ReadIn(Input16, Buffer16);
			IOHandler32.ReadIn(Input32, Buffer32);

			ArrayXXf BufferTest;
			BufferTest.resize(BufferSize, Channels);
			for (auto i = 0; i < Channels; i++)
			{
				BufferTest.col(i) = ArrayXf::LinSpaced(BufferSize, (float)(i) / (Channels * BufferSize), (float)(i + (BufferSize - 1) * Channels) / (Channels * BufferSize));
			}
			float error = (Buffer16 - BufferTest).abs2().sum();
			Assert::IsTrue(error < 1e-10f);
			error = (Buffer32 - BufferTest).abs2().sum();
			Assert::IsTrue(error < 1e-10f);
			outputLog << "Test succesful." << std::endl;
		}

		TEST_METHOD(ReadWriteVoidPtr)
		{
			outputLog << "Running FormatHandlingTest->ReadWriteVoidPtr.\n";
			FormatHandling<float> IOHandler;
			const int Channels = 5;
			const int BufferSize = 100;
			ArrayXXf Buffer;
			Buffer.resize(BufferSize, Channels);
			float Input[Channels * BufferSize];
			for (auto i = 0; i < Channels * BufferSize; i++) { Input[i] = static_cast<float>(i); }

			void* InputVoid = Input; // assume input is coming from some C-code as a void pointer
			IOHandler.ReadIn(InputVoid, Buffer);

			ArrayXXf BufferTest;
			BufferTest.resize(BufferSize, Channels);
			for (auto i = 0; i < Channels; i++)
			{
				BufferTest.col(i) = ArrayXf::LinSpaced(BufferSize, (float)i, (float)(i + (BufferSize - 1) * Channels));
			}
			float error = (Buffer - BufferTest).abs2().sum();
			Assert::IsTrue(error < 1e-10f);

			float Output[Channels * BufferSize];
			void* OutputVoid = Output; // assume output is given as a void ptr
			IOHandler.WriteOut(Buffer, OutputVoid);
			error = 0.f;
			for (auto i = 0; i < Channels * BufferSize; i++) { error += fabs(Output[i] - Input[i]); }
			Assert::IsTrue(error < 1e-10f);
			outputLog << "Test succesful." << std::endl;
		}
		TEST_METHOD(BaseSmartPointer)
		{
			outputLog << "Running FormatHandlingTest->BaseSmartPointer.\n";
			std::unique_ptr<FormatHandlingBase> IOHandler = std::make_unique<FormatHandling<float>>();
			float input[500];
			for (auto i = 0; i < 500; i++) { input[i] = static_cast<float>(i); }
			ArrayXXf buffer(100, 5);
			auto start = std::chrono::steady_clock::now();
			for (int i = 0; i < 100; i++)
			{
				IOHandler->ReadIn(input, buffer);
			}
			auto end = std::chrono::steady_clock::now();
			outputLog << "Elapsed Time base pointer: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << " nanoseconds\n";

			FormatHandling<float> IOHandler2;
			ArrayXXf buffer2(100, 5);
			start = std::chrono::steady_clock::now();
			for (int i = 0; i < 100; i++)
			{
				IOHandler2.ReadIn(input, buffer2);
			}
			end = std::chrono::steady_clock::now();
			outputLog << "Elapsed Time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << " nanoseconds\n";

			float error = (buffer - buffer2).abs2().sum();
			Assert::IsTrue(error < 1e-10f);
			outputLog << "Test succesful." << std::endl;
		}
	};
}