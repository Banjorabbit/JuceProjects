#include "unittest.h"
#include "../src/AsynchronousBase/Distortion.h"
#include "../src/AsynchronousBase/FFTConvolution.h"
#include "../src/AsynchronousBase/IIR2ndDF.h"
#include "../src/AsynchronousBase/IIR2ndNonLin.h"
#include "../src/AsynchronousBase/IIR2ndTimeVarying.h"
#include "../src/AsynchronousBase/LimiterHard.h"
#include "../src/AsynchronousBase/NonparametricEqualizer.h"
#include "../src/AsynchronousBase/PipelinePreProcessing.h"
#include "../src/AsynchronousBase/PitchShift.h"
#include "../src/AsynchronousBase/PitchShiftAdaptiveResolution.h"
#include "../src/AsynchronousBase/PitchShiftPhaseLocking.h"
#include "../src/AsynchronousBase/SeparateTonal.h"
#include "../src/AsynchronousBase/SeparateTonalTextureTransient.h"
#include "../src/AsynchronousBase/SeparateTransient.h"
#include "../src/AsynchronousBase/VirtualizationHeadphones.h"
#include "../Utilities/AudioFile.h"
#include "../Utilities/FormatHandling.h"
#include "../Utilities/pffft.h"

using namespace InterfaceTests;


namespace AsynchronousBaseTests
{
	LoggerOutput outputLog;

	TEST_CLASS(DistortionTest)
	{
		TEST_METHOD(InterfaceAsynchronous)
		{
			outputLog << "Running DistortionTest->InterfaceAsynchronous.\n";
			auto flag = AlgorithmInterfaceAsynchronousTest<Distortion5thOrderOdd>();
			Assert::IsTrue(flag);
		}
	};

	TEST_CLASS(NonparametricEqualizerTest)
	{
		TEST_METHOD(InterfaceAsynchronous)
		{
			outputLog << "Running NonparametricEqualizerTest->InterfaceAsynchronous.\n";
			ArrayXf frequencies(5);
			ArrayXf gaindB(5);
			frequencies << 80, 500, 1000, 2400, 4000;
			gaindB << 5, 10, -10, 0, 7;
			I::NonparametricEqualizerPersistent inputX = { frequencies, gaindB };
			auto flag = AlgorithmInterfaceAsynchronousTest<NonparametricEqualizer>(inputX);
			Assert::IsTrue(flag);
		}
	};

	TEST_CLASS(PipelinePreProcessingTest)
	{
		TEST_METHOD(InterfaceAsynchronous)
		{
			outputLog << "Running PipelinePreProcessingTest->Interface.\n";
			auto flag = AlgorithmInterfaceAsynchronousTest<PipelinePreProcessing>();
			Assert::IsTrue(flag);
		}
	};

	TEST_CLASS(PitchShiftTest)
	{
		TEST_METHOD(InterfaceAsynchronous)
		{
			auto flag = AlgorithmInterfaceAsynchronousTest<PitchShift>();
			Assert::IsTrue(flag);
		}

	};

	TEST_CLASS(PitchShiftAdaptiveResolutionTest)
	{
		TEST_METHOD(InterfaceAsynchronous)
		{
			auto flag = AlgorithmInterfaceAsynchronousTest<PitchShiftAdaptiveResolution>();
			Assert::IsTrue(flag);
		}

	};

	TEST_CLASS(PitchShiftPhaseLockingTest)
	{
		TEST_METHOD(InterfaceAsynchronous)
		{
			auto flag = AlgorithmInterfaceAsynchronousTest<PitchShiftPhaseLocking>();
			Assert::IsTrue(flag);
		}

	};

	TEST_CLASS(SeparateTonalTest)
	{
		TEST_METHOD(InterfaceAsynchronous)
		{
			auto flag = AlgorithmInterfaceAsynchronousTest<SeparateTonal>();
			Assert::IsTrue(flag);
		}

		TEST_METHOD(LatencyTest)
		{
			SeparateTonal tonalSeparator;
			auto c = tonalSeparator.GetCoefficients();
			c.NChannelsIn = 1;
			tonalSeparator.Initialize(c);
			auto latency = tonalSeparator.GetLatencySamples();
			
			ArrayXf input(c.BufferSize * 20), output(c.BufferSize * 20);
			input.setRandom();
			for (auto i = 0; i < 20; i++)
			{
				tonalSeparator.Process(input.segment(i*c.BufferSize, c.BufferSize), output.segment(i*c.BufferSize, c.BufferSize));
			}
			
			outputLog << "Latency: " << latency << "\n";
			outputLog << "Diff: " << (output.tail(20*c.BufferSize - latency) - input.head(20*c.BufferSize - latency)).abs2().sum() / input.head(20*c.BufferSize - latency).abs2().sum() << std::endl;
		}
	};

	TEST_CLASS(SeparateTonalTextureTransientTest)
	{
		TEST_METHOD(InterfaceAsynchronous)
		{
			auto flag = AlgorithmInterfaceAsynchronousTest<SeparateTonalTextureTransient>();
			Assert::IsTrue(flag);
		}
	};

	TEST_CLASS(SeparateTransientTest)
	{
		TEST_METHOD(InterfaceAsynchronous)
		{
			auto flag = AlgorithmInterfaceAsynchronousTest<SeparateTransient>();
			Assert::IsTrue(flag);
		}
	};

	TEST_CLASS(IIR2ndDFTest)
	{
		TEST_METHOD(InterfaceAsynchronous2ndCascaded)
		{
			IIR2ndCascaded filter;
			auto c = filter.GetCoefficients();
			ArrayXXf sos(c.Nsos, 6);
			sos.setRandom();
			AlgorithmInterfaceAsynchronousTest<IIR2ndCascaded>(sos);
		}

		TEST_METHOD(CheckCalculationIIR2ndCascaded)
		{
			IIR2ndCascaded filter;
			int bufferSize = 16;
			auto c = filter.GetCoefficients();
			c.NChannelsIn = 1;
			c.Nsos = 3;
			c.BufferSize = bufferSize;
			filter.Initialize(c);
			ArrayXXf sos(c.Nsos, 6);
			sos << 1.f, 1.f, -1.f, 1.f, .2f, -.2f, 1.f, -1.f, -1.f, 1.f, .3f, -.4f, 1.f, 2.f, 1.f, 1.f, -.1f, .05f;
			filter.SetPersistentInput(sos);
			ArrayXf input(bufferSize), output(bufferSize), outputRef(bufferSize);
			input.setZero();
			input(0) = 1.f;
			outputRef << 1.000000000000000f, 1.600000000000000f, -2.100000000000000f, -4.234999999999999f, -1.409400000000000f, 0.056560000000000f, 0.118665000000000f, 0.142948000000000f, 0.049654860000000f, 0.056578021000000f, 0.011437689000000f, 0.020982221200000f, -0.000422592759000f, 0.008579439073600f, -0.002496334514100f, 0.004145097177685f;
			filter.Process(input, output);
			float error = (output - outputRef).abs2().sum() / outputRef.abs2().sum();
			outputLog << "error: " << error << std::endl;
			Assert::IsTrue(error < 1e-10f);
		}
		TEST_METHOD(InterfaceAsynchronous2ndOrder)
		{
			outputLog << "Running IIRFiltersTest->InterfaceAsynchronous2ndOrder.\n";
			auto flag = AlgorithmInterfaceAsynchronousTest<IIR2ndDF2>();
			Assert::IsTrue(flag);

			flag = AlgorithmInterfaceAsynchronousTest<IIR2ndDF2Transposed>();
			Assert::IsTrue(flag);
		}

		TEST_METHOD(InterfaceAsynchronous2ndOrderComplex)
		{
			outputLog << "Running IIRFiltersTest->InterfaceAsynchronous2ndDF2Complex.\n";
			auto flag = AlgorithmInterfaceAsynchronousTest<IIR2ndDF2Complex>();
			Assert::IsTrue(flag);

			flag = AlgorithmInterfaceAsynchronousTest<IIR2ndDF2TransposedComplex>();
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

			ArrayXXf input(bufferSize, c.NChannelsIn), outputI(bufferSize, c.NChannelsIn), outputII(bufferSize, c.NChannelsIn);
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
		TEST_METHOD(AsynchronousBaseTest)
		{
			IIR2ndNonLinLowpass filter;
			outputLog << "Un-initialized.\n";
			outputLog << "GetInitializedAsynchronous(): " << filter.GetInitializedAsynchronous() << "\n";
			outputLog << "GetBufferSizeExternal(): " << filter.GetBufferSizeExternal() << "\n";
			outputLog << "GetSynchronousProcessing(): " << filter.GetSynchronousProcessing() << "\n";
			outputLog << "GetNChannelsIn(): " << filter.GetNChannelsIn() << "\n";
			outputLog << "GetNChannelsOut(): " << filter.GetNChannelsOut() << "\n";
			outputLog << "GetBufferSize(): " << filter.GetBufferSize() << "\n";
			outputLog << "GetLatencySamples(): " << filter.GetLatencySamples() << "\n";

			filter.Initialize();
			outputLog << "Initialized.\n";
			outputLog << "GetInitializedAsynchronous(): " << filter.GetInitializedAsynchronous() << "\n";
			outputLog << "GetBufferSizeExternal(): " << filter.GetBufferSizeExternal() << "\n";
			outputLog << "GetSynchronousProcessing(): " << filter.GetSynchronousProcessing() << "\n";
			outputLog << "GetNChannelsIn(): " << filter.GetNChannelsIn() << "\n";
			outputLog << "GetNChannelsOut(): " << filter.GetNChannelsOut() << "\n";
			outputLog << "GetBufferSize(): " << filter.GetBufferSize() << "\n";
			outputLog << "GetLatencySamples(): " << filter.GetLatencySamples() << "\n";

			filter.InitializeAsynchronous(11);
			outputLog << "Initialized Asynchronous.\n";
			outputLog << "GetInitializedAsynchronous(): " << filter.GetInitializedAsynchronous() << "\n";
			outputLog << "GetBufferSizeExternal(): " << filter.GetBufferSizeExternal() << "\n";
			outputLog << "GetSynchronousProcessing(): " << filter.GetSynchronousProcessing() << "\n";
			outputLog << "GetNChannelsIn(): " << filter.GetNChannelsIn() << "\n";
			outputLog << "GetNChannelsOut(): " << filter.GetNChannelsOut() << "\n";
			outputLog << "GetBufferSize(): " << filter.GetBufferSize() << "\n";
			outputLog << "GetLatencySamples(): " << filter.GetLatencySamples() << "\n";

			auto c = filter.GetCoefficients();
			c.BufferSize = 256;
			filter.InitializeAsynchronous(c, 128);
			outputLog << "Initialized Asynchronous with BufferSize=256.\n";
			outputLog << "GetInitializedAsynchronous(): " << filter.GetInitializedAsynchronous() << "\n";
			outputLog << "GetBufferSizeExternal(): " << filter.GetBufferSizeExternal() << "\n";
			outputLog << "GetSynchronousProcessing(): " << filter.GetSynchronousProcessing() << "\n";
			outputLog << "GetNChannelsIn(): " << filter.GetNChannelsIn() << "\n";
			outputLog << "GetNChannelsOut(): " << filter.GetNChannelsOut() << "\n";
			outputLog << "GetBufferSize(): " << filter.GetBufferSize() << "\n";
			outputLog << "GetLatencySamples(): " << filter.GetLatencySamples() << "\n";

			outputLog << std::endl;
		}
		TEST_METHOD(InterfaceAsynchronousNonLinLowpass)
		{
			outputLog << "Running IIRFiltersTest->InterfaceAsynchronousNonLinLowpass." << std::endl;
			auto flag = AlgorithmInterfaceAsynchronousTest<IIR2ndNonLinLowpass>();
			Assert::IsTrue(flag);
		}

		TEST_METHOD(InterfaceAsynchronousNonLinBandpass)
		{
			outputLog << "Running IIRFiltersTest->InterfaceAsynchronousNonLinBandpass.\n";
			auto flag = AlgorithmInterfaceAsynchronousTest<IIR2ndNonLinBandpass>();
			Assert::IsTrue(flag);
		}

		TEST_METHOD(CalculationNonLinLowpass)
		{
			IIR2ndNonLinLowpass filter;
			auto s = filter.GetSetup();
			int bufferSize = 5;
			s.Coefficients.SampleRate = 44100.f;
			s.Coefficients.NChannelsIn = 1;
			s.Coefficients.BufferSize = bufferSize;
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
		TEST_METHOD(InterfaceAsynchronousTimeVarying)
		{
			outputLog << "Running IIRFiltersTest->InterfaceAsynchronousTimeVarying.\n";
			auto flag = AlgorithmInterfaceAsynchronousTest<IIR2ndTimeVaryingLowPassFilter>();
			Assert::IsTrue(flag);

			flag = AlgorithmInterfaceAsynchronousTest<IIR2ndTimeVaryingHighPassFilter>();
			Assert::IsTrue(flag);

			flag = AlgorithmInterfaceAsynchronousTest<IIR2ndTimeVaryingBandPassFilter>();
			Assert::IsTrue(flag);

			flag = AlgorithmInterfaceAsynchronousTest<IIR2ndTimeVaryingNotchFilter>();
			Assert::IsTrue(flag);

			flag = AlgorithmInterfaceAsynchronousTest<IIR2ndTimeVaryingBellFilter>();
			Assert::IsTrue(flag);

			flag = AlgorithmInterfaceAsynchronousTest<IIR2ndTimeVaryingHighShelfFilter>();
			Assert::IsTrue(flag);

			flag = AlgorithmInterfaceAsynchronousTest<IIR2ndTimeVaryingLowShelfFilter>();
			Assert::IsTrue(flag);
		}
	};

	TEST_CLASS(LimiterHardTest)
	{
		TEST_METHOD(InterfaceAsynchronous)
		{
			outputLog << "Running LimiterHardTest->InterfaceAsynchronous.\n";
			auto flag = AlgorithmInterfaceAsynchronousTest<LimiterHard>();
			Assert::IsTrue(flag);
		}
		// delay input signal and save result to file
		TEST_METHOD(DelaylinearArray)
		{
			LimiterHard limiter;
			auto c = limiter.GetCoefficients();
			c.BufferSize = 64;
			c.HoldTimeMS = 7.f;
			c.NChannelsIn = 2;
			c.SampleRate = 44100.f;
			limiter.Initialize(c);

			AudioFile<float> audioFile, audioFileRef;
			audioFile.setNumChannels(c.NChannelsIn);
			audioFile.setBitDepth(24);
			audioFile.setSampleRate(static_cast<int>(c.SampleRate));
			audioFileRef.setNumChannels(c.NChannelsIn);
			audioFileRef.setBitDepth(24);
			audioFileRef.setSampleRate(static_cast<int>(c.SampleRate));

			int nBuffers = static_cast<int>(c.SampleRate/c.BufferSize*c.HoldTimeMS*5/1000);
			ArrayXXf input(nBuffers * c.BufferSize, c.NChannelsIn);
			ArrayXXf output(c.BufferSize, c.NChannelsIn);
			input.col(0).setLinSpaced(nBuffers*c.BufferSize, 0, 1);
			input.col(1).setLinSpaced(nBuffers*c.BufferSize, 0, 1);
			for (auto i = 0; i < nBuffers; i++)
			{
				limiter.Process(input.middleRows(i*c.BufferSize, c.BufferSize), output);
				for (auto j = 0; j < c.BufferSize; j++)
				{
					for (auto channel = 0; channel < c.NChannelsIn; channel++)
					{
						audioFile.samples[channel].push_back(output(j, channel));
						audioFileRef.samples[channel].push_back(input.middleRows(i*c.BufferSize, c.BufferSize)(j, channel));
					}
				}
			}
			outputLog << "Saving file.\n";
			audioFile.save("../../../../Temp/LimiterHard.wav", AudioFileFormat::Wave); // truncate to 24bits 
			audioFileRef.save("../../../../Temp/LimiterHardRef.wav", AudioFileFormat::Wave); // truncate to 24bits 
		}
	};

	TEST_CLASS(FFTConvolutionTest)
	{
		TEST_METHOD(InterfaceAsynchronous)
		{
			FFTConvolution conv;
			auto c = conv.GetCoefficients();
			ArrayXf filter(c.FilterSize);
			filter.setRandom();
			auto flag = AlgorithmInterfaceAsynchronousTest<FFTConvolution>(filter);
			Assert::IsTrue(flag);
		}
	};

	TEST_CLASS(VirtualizationHeadphonesTest)
	{
		TEST_METHOD(InterfaceAsynchronous)
		{
			outputLog << "Running VirtualizationHeadphonesTest->Interface.\n";
			auto flag = AlgorithmInterfaceAsynchronousTest<VirtualizationHeadphones>();
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
			ArrayXXf inputx1(c.BufferSize, c.NChannelsIn), output(c.BufferSize, 2);
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
			ArrayXXf inputx2(cx2.BufferSize, cx2.NChannelsIn), outputx2(cx2.BufferSize, 2);
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
			ArrayXXf inputx4(cx4.BufferSize, cx4.NChannelsIn), outputx4(cx4.BufferSize, 2);
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

			ArrayXXf input(c.BufferSize, c.NChannelsIn), output(c.BufferSize, 2);
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