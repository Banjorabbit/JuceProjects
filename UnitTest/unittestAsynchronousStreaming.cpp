#include "unittest.h"
#include "../src/AsynchronousStreaming/PipelinePreProcessing.h"
#include "../src/AsynchronousStreaming/FFTConvolution.h"
#include "../src/AsynchronousStreaming/IIR2ndDF.h"
#include "../src/AsynchronousStreaming/IIR2ndNonLin.h"
#include "../src/AsynchronousStreaming/IIR2ndTimeVarying.h"
#include "../src/AsynchronousStreaming/NonparametricEqualizer.h"
#include "../src/AsynchronousStreaming/PitchShift.h"
#include "../src/AsynchronousStreaming/PitchShiftAdaptiveResolution.h"
#include "../src/AsynchronousStreaming/VirtualizationHeadphones.h"
#include "../Utilities/AudioFile.h"
#include "../Utilities/FormatHandling.h"
#include "../Utilities/pffft.h"

using namespace InterfaceTests;


namespace AsynchronousStreamingTests
{
	LoggerOutput outputLog;

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
			I::NonparametricEqualizerPersistent inputX = { frequencies, gaindB };
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

	TEST_CLASS(PitchShiftAdaptiveResolutionTest)
	{
		TEST_METHOD(InterfaceStreaming)
		{
			PitchShiftAdaptiveResolution timeStretch;
			auto c = timeStretch.GetCoefficients();
			ArrayXXf input(c.BufferSize, 1);
			input.setRandom();
			ArrayXXf output(c.BufferSize, 1);
			auto sampleRate = 48e3f;
			auto flag = AlgorithmInterfaceStreamingTest<PitchShiftAdaptiveResolutionStreaming>(input, output, sampleRate);
			Assert::IsTrue(flag);
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
				error += std::abs(outputx2(2 * i, 0) * 2 - output(i, 0));
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
				error += std::abs(outputx4(4 * i, 0) * 4 - output(i, 0));
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

			ArrayXXf input(c.BufferSize, c.NChannels), output(c.BufferSize, 2);
			input.setZero();
			virtualHeadphones.Process(input, output);
			float error = (output.bottomRows(c.BufferSize - delay) - input.rightCols<1>().head(c.BufferSize - delay).replicate(1, 2)).abs2().sum();
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
}