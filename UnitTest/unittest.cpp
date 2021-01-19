#include "unittest.h"
#include "../src/CircBuffer.h"
#include "../src/CubicSpline.h"
#include "../src/DesignIIRMinPhase.h"
#include "../src/DesignIIRNonParametric.h"
#include "../src/DesignFIRNonParametric.h"
#include "../src/FFT.h"
#include "../src/Filterbank.h"
#include "../src/FilterMinMax.h"
#include "../src/FIRMinPhase.h"
#include "../src/InterpolationCubic.h"
#include "../src/SpectralMass.h"
#include "../src/StateVariableFilter.h"
#include "../src/ToeplitzSolver.h"
#include "../src/Upsampling2XCubic.h"
#include "../Utilities/AudioFile.h"

using namespace InterfaceTests;

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

namespace Algorithms
{		
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

	TEST_CLASS(FilterMinMaxTest)
	{
		TEST_METHOD(InterfaceStreamingMinMax)
		{
			StreamingMinMax streaming;
			auto c = streaming.GetCoefficients();
			ArrayXXf input(1000, c.NChannels);
			input.setRandom();
			ArrayXXf MinValues(1000, c.NChannels);
			ArrayXXf MaxValues(1000, c.NChannels);
			O::FilterMinMax output = { MinValues, MaxValues };
			auto flag = AlgorithmInterfaceTest<StreamingMinMax>(input, output);
			Assert::IsTrue(flag);
		}

		TEST_METHOD(InterfaceStreamingMax)
		{
			StreamingMax streaming;
			auto c = streaming.GetCoefficients();
			ArrayXXf input(1000, c.NChannels);
			input.setRandom();
			ArrayXXf output(1000, c.NChannels);
			auto flag = AlgorithmInterfaceTest<StreamingMax>(input, output);
			Assert::IsTrue(flag);
		}

		TEST_METHOD(InterfaceFilterMinMax)
		{
			FilterMinMax filter;
			auto c = filter.GetCoefficients();
			ArrayXXf input(1000, c.NChannels);
			input.setRandom();
			ArrayXXf MinValues(1000, c.NChannels);
			ArrayXXf MaxValues(1000, c.NChannels);
			O::FilterMinMax output = { MinValues, MaxValues };
			auto flag = AlgorithmInterfaceTest<FilterMinMax>(input, output);
			Assert::IsTrue(flag);

			flag = AlgorithmInterfaceTest<FilterMax>(input, MaxValues);
			Assert::IsTrue(flag);

			flag = AlgorithmInterfaceTest<FilterMin>(input, MinValues);
			Assert::IsTrue(flag);
		}

		TEST_METHOD(CheckCalculationMaxMin)
		{
			// streaming algorithms
			StreamingMinMax streaming;
			auto c = streaming.GetCoefficients();
			c.NChannels = 1;
			c.Length = 5;
			streaming.Initialize(c);

			StreamingMax sMax;
			auto cMax = sMax.GetCoefficients();
			cMax.NChannels = 1;
			cMax.Length = 5;
			sMax.Initialize(cMax);

			StreamingMin sMin;
			auto cMin = sMin.GetCoefficients();
			cMin.NChannels = 1;
			cMin.Length = 5;
			sMin.Initialize(cMin);

			ArrayXf input(20);
			input.setRandom();
			ArrayXf minValue(20), maxValue(20), output(20);
			streaming.Process(input, { minValue, maxValue });
			sMax.Process(input, output);

			float error = (output - maxValue).abs2().sum();
			outputLog << "Error: " << error << std::endl;
			Assert::IsTrue(error == 0);

			sMin.Process(input, output);
			error = (output - minValue).abs2().sum();
			outputLog << "Error: " << error << std::endl;
			Assert::IsTrue(error == 0);

			// filter algorithms
			FilterMinMax filter;
			auto cf = filter.GetCoefficients();
			cf.NChannels = 1;
			cf.Length = 5;
			filter.Initialize(cf);

			FilterMax fMax;
			auto cfMax = fMax.GetCoefficients();
			cfMax.NChannels = 1;
			cfMax.Length = 5;
			fMax.Initialize(cfMax);

			FilterMin fMin;
			auto cfMin = fMin.GetCoefficients();
			cfMin.NChannels = 1;
			cfMin.Length = 5;
			fMin.Initialize(cfMin);

			filter.Process(input, { minValue, maxValue });
			fMax.Process(input, output);

			error = (output - maxValue).abs2().sum();
			outputLog << "Error: " << error << std::endl;
			Assert::IsTrue(error == 0);

			fMin.Process(input, output);
			error = (output - minValue).abs2().sum();
			outputLog << "Error: " << error << std::endl;
			Assert::IsTrue(error == 0);
		}

		TEST_METHOD(CheckCalculation)
		{
			StreamingMinMax streaming;
			auto c = streaming.GetCoefficients();
			c.NChannels = 1;
			c.Length = 5;
			streaming.Initialize(c);

			ArrayXf input(20);
			input << 1.f, 1.2f, 1.5f, 1.8f, 1.4f, 1.15f, 0.6f, 2.1f, 1.3f, 1.05f, 1.9f, 2.2f, 1.3f, 0.7f, 0.3f, 1.4f, 1.7f, 1.8f, 1.41f, 1.2f;
			ArrayXf minValue(20), maxValue(20), minRef(20), maxRef(20);
			minRef << 1.f, 1.f, 1.f, 1.f, 1.f, 1.15f, 0.6f, 0.6f, 0.6f, 0.6f, 0.6f, 1.05f, 1.05f, 0.7f, 0.3f, 0.3f, 0.3f, 0.3f, 0.3f, 1.2f;
			maxRef << 1.f, 1.2f, 1.5f, 1.8f, 1.8f, 1.8f, 1.8f, 2.1f, 2.1f, 2.1f, 2.1f, 2.2f, 2.2f, 2.2f, 2.2f, 2.2f, 1.7f, 1.8f, 1.8f, 1.8f;
			O::FilterMinMax output = { minValue, maxValue };
			
			// use ResetInitializeValue to set start memory
			streaming.ResetInitialValue(input(0));
			streaming.Process(input, output);
			float error = (output.MaxValue - maxRef).abs2().sum() + (output.MinValue - minRef).abs2().sum();
			outputLog << "Error: " << error << std::endl;
			Assert::IsTrue(error == 0);

			// repeat with FilterMinMax and check it gives the same values shifted
			FilterMinMax filter;
			auto cF = filter.GetCoefficients();
			cF.NChannels = 1;
			cF.Length = 5;
			filter.Initialize(cF);
			ArrayXf minValue2(20), maxValue2(20);
			O::FilterMinMax output2 = { minValue2, maxValue2 };

			filter.Process(input, output2);

			error = (output.MaxValue.bottomRows(18) - output2.MaxValue.topRows(18)).abs2().sum() + (output.MinValue.bottomRows(18) - output2.MinValue.topRows(18)).abs2().sum();
			outputLog << "error: " << error << std::endl;
			Assert::IsTrue(error == 0);
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

	TEST_CLASS(SpectralMassTest)
	{
		TEST_METHOD(Interface)
		{
			outputLog << "Running SpectralMassTest->Interface.\n";
			SpectralMass spectralMass;
			auto c = spectralMass.GetCoefficients();
			ArrayXXf input(128, c.NChannels);
			input.setRandom();
			ArrayXXf output(128, c.NChannels);
			auto flag = AlgorithmInterfaceTest<SpectralMass>(input, output);
			Assert::IsTrue(flag);
		}

		TEST_METHOD(CalculationTest)
		{
			outputLog << "Running SpectralMassTest->CalculationTest.\n";
			SpectralMass spectralMass;
			auto s = spectralMass.GetSetup();
			s.Coefficients.NChannels = 1;
			s.Coefficients.SampleRate = 44100.f;
			s.Parameters.SmoothTConstant = 1e-40f;
			spectralMass.Initialize(s);
			ArrayXf input(16), output(16), outputRef(16);
			input << -1.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, 1.f, 1.f, -1.f, -1.f, 1.f, 1.f, -1.f, -1.f, 1.f;
			outputRef << 1000, 22050, 22050, 22050, 22050, 22050, 22050, 22050, 22050, 11025, 11025, 11025, 11025, 11025, 11025, 11025;
			spectralMass.Process(input, output);
			float error = (output - outputRef).abs2().sum();
			Assert::IsTrue(error == 0.f);

			// run again
			spectralMass.Process(input, output);
			outputRef(0) = 22050;
			error = (output - outputRef).abs2().sum();
			Assert::IsTrue(error == 0.f);
		}
	};

	TEST_CLASS(StateVariableFilterTest)
	{
		TEST_METHOD(Interface)
		{
			StateVariableFilter filter;
			auto c = filter.GetCoefficients();
			auto bufferSize = 128;
			ArrayXXf input(bufferSize, c.NChannels), BandPass(bufferSize, c.NChannels), HighPass(bufferSize, c.NChannels), LowPass(bufferSize, c.NChannels);
			input.setRandom();
			O::StateVariableFilter output = { HighPass, LowPass, BandPass };
			auto flag = AlgorithmInterfaceTest<StateVariableFilter>(input, output);
			Assert::IsTrue(flag);
		}
	};

	TEST_CLASS(ToeplitzSolverTest)
	{
		TEST_METHOD(Interface)
		{
			outputLog << "Running ToeplitzSolverTest->Interface.\n";
			ToeplitzSolver toeplitzSolver;
			ArrayXcf aToeplitz(32);
			ArrayXXcf bRightHand(32, 8);
			I::ToeplitzSolver input = { aToeplitz, bRightHand };
			ArrayXXcf output(32, 8);
			auto flag = AlgorithmInterfaceTest<ToeplitzSolver>(input, output);
			Assert::IsTrue(flag);
		}

		TEST_METHOD(CalculationTest)
		{
			using namespace std::literals;

			outputLog << "running ToeplitzSolverTest->CalculationTest.\n";
			ToeplitzSolver toeplitzSolver;
			ArrayXcf aToeplitz(2);
			ArrayXXcf bRightHand(2, 3);
			aToeplitz << 0.4882f - 0.1961if, -0.1774f + 1.4193if;
			bRightHand << -0.2103f + 0.2486if, 0.6536f - 0.0990if, -0.7586f + 0.6976if, -0.3118f - 0.2045if, -0.7980f - 0.8078if, -0.3473f - 0.6120if;
			I::ToeplitzSolver input = { aToeplitz, bRightHand };
			ArrayXXcf output(2, 3);
			ArrayXXcf outputRef(2, 3);
			outputRef << 0.1825f - 0.3271if, 0.4846f - 0.4906if, 0.5752f - 0.5329if, 0.3284f + 0.1247if, 0.1188f - 0.3764if, 0.8237f + 0.5558if;
			toeplitzSolver.Initialize();
			toeplitzSolver.Process(input, output);
			float error = (output - outputRef).abs2().sum();
			outputLog << "Error: " << error << std::endl;
			Assert::IsTrue(error < 1e-7f);
		}
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

		TEST_METHOD(UpsampleImpulse)
		{
			outputLog << "Running Upsampling2XCubicTest->UpsampleImpulse.\n";
			UpsamplingPower2Cubic upsampling;
			auto p = upsampling.GetParameters();
			p.UpsamplingFactorPow2 = 3;
			upsampling.SetParameters(p);
			upsampling.Initialize();

			ArrayXf input(8);
			input.setZero();
			input(3) = 1;
			ArrayXf output(8 * 8);
			upsampling.Process(input, output);
			outputLog << "Output: " << output << std::endl;
		}
	};
}