#include "unittest.h"
#include "../src/FrequencyDomain/BeamformerAdaptive.h"
#include "../src/FrequencyDomain/CriticalBands.h"
#include "../src/FrequencyDomain/DetectTonal.h"
#include "../src/FrequencyDomain/DetectTransient.h"
#include "../src/FrequencyDomain/DetectVoiceActivation.h"
#include "../src/FrequencyDomain/EchoCancellerMomentum.h"
#include "../src/FrequencyDomain/EchoSuppressionCovariance.h"
#include "../src/FrequencyDomain/GainCalculation.h"
#include "../src/FrequencyDomain/InterpolationTemporal.h"
#include "../src/FrequencyDomain/LoudnessModel.h"
#include "../src/FrequencyDomain/MinPhaseSpectrum.h"
#include "../src/FrequencyDomain/NoiseEstimation.h"
#include "../src/FrequencyDomain/NoiseSuppression.h"
#include "../Utilities/AudioFile.h"

using namespace InterfaceTests;

namespace FrequencyDomainTests
{
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
				Assert::IsTrue(std::abs(Data->Rx[i](1, 0)) > 10 * std::abs(Data->Rn[i](1, 0)));
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
			DetectVoiceActivation VAD;
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
			DetectVoiceActivation voiceActivationDetection;
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
			ArrayXXf input(100, 10);
			input.setRandom();
			ArrayXXf output(100, 10);
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
			ArrayX2f xEnergy(c.NBands, 2);
			ArrayX2f xPhase(c.NBands, 2);
			float point = 0.2f;
			xEnergy.setRandom().abs();
			xPhase.setRandom();
			I::InterpolationTemporal input = { xEnergy, xPhase, point };
			ArrayXcf output(c.NBands);
			auto flag = AlgorithmInterfaceTest<InterpolationTemporal>(input, output);
			Assert::IsTrue(flag);

		}
	};

	TEST_CLASS(LoudnessModelTest)
	{
		TEST_METHOD(Interface)
		{
			LoudnessModel model;
			model.Initialize();
			auto NBandsCritical = model.ConvertBands.GetNBandsCritical();
			auto c = model.GetCoefficients();
			ArrayXXf input(c.NBands, c.NChannels);
			input.setRandom();
			input = input.abs();
			ArrayXXf output(NBandsCritical, c.NChannels);
			auto flag = AlgorithmInterfaceTest<LoudnessModel>(input, output);
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

	TEST_CLASS(DetectTonalTest)
	{
		TEST_METHOD(Interface)
		{
			DetectTonal tonalDetector;
			auto c = tonalDetector.GetCoefficients();
			ArrayXXcf input(c.NBands, c.NChannels);
			input.setRandom();
			input *= 1;
			Array<bool, Dynamic,Dynamic> output(c.NBands, c.NChannels);
			auto flag = AlgorithmInterfaceTest<DetectTonal>(input, output);
			Assert::IsTrue(flag);
		}
		
		// Algorithm changed to output bools, so this test is no longer supported
		//TEST_METHOD(CheckCalculation)
		//{
		//	DetectTonal tonalDetector;
		//	auto setup = tonalDetector.GetSetup();
		//	setup.Coefficients.FilterbankRate = 5512.5f;
		//	setup.Coefficients.FilterRelativeRange = 0.1875f;
		//	setup.Coefficients.NBands = 17;
		//	setup.Parameters.TonalTC = 0.02f;
		//	setup.Parameters.TonalTransientRatio = 0.262853538008622f;
		//	tonalDetector.Initialize(setup);
		//	ArrayXf input(17);
		//	input << 6.101540931246211f, 2.157399825081734f, 0.030358307769100f, 0.001622377620863f, 0.000002650702480f, 0.000082119432769f, 0.000023003709579f, 0.000142849128474f, 0.000095319771848f, 0.000007427525950f, 0.000011033106498f, 0.000002778327015f, 0.000001644832597f, 0.000002744397722f, 0.000003318713830f, 0.000000526210865f, 0.000000215532469f;
		//	ArrayXf output(17), ref(17);
		//	tonalDetector.Process(input, output);
		//	ref << 0.003980699139697f, 0.003660735591129f, 0.001949118358953f, 0.000054814734024f, 0.000014509997412f, 0.000076508332610f, 0.001308871015829f, 0.002831560773398f, 0.002603381025513f, 0.002606893096354f, 0.000543286755640f, 0.001059095995554f, 0.000708450435246f, 0.000726312436408f, 0.002118427013139f, 0.002427550790790f, 0.000517556810844f;
		//	float error = (ref - output).abs2().sum() / ref.abs2().sum();
		//	Assert::IsTrue(error < 1e-10f);
		//}
	};

	TEST_CLASS(DetectTransientTest)
	{
		TEST_METHOD(Interface)
		{
			outputLog << "Running DetectTransientTest->Interface.\n";
			DetectTransient TransientDetection;
			TransientDetection.Initialize();
			auto c = TransientDetection.GetCoefficients();
			ArrayXXf input(c.NBands, c.NChannels);
			int nBandsCritical = TransientDetection.ConvertBands.GetNBandsCritical();
			Array<bool, Dynamic, Dynamic> output(nBandsCritical, c.NChannels);
			auto flag = AlgorithmInterfaceTest<DetectTransient>(input, output);
			Assert::IsTrue(flag);
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

	TEST_CLASS(CriticalBandsTest)
	{
		TEST_METHOD(Interface)
		{
			CriticalBands model;
			model.Initialize();
			auto NBandsCritical = model.GetNBandsCritical();
			auto c = model.GetCoefficients();
			ArrayXXf input(c.NBands, 2);
			input.setRandom();
			input = input.abs();
			ArrayXXf output(NBandsCritical, 2);
			auto flag = AlgorithmInterfaceTest<CriticalBands>(input, output);
			Assert::IsTrue(flag);
		}

		TEST_METHOD(GetCenterFrequenciesTest)
		{
			CriticalBands model;
			int Nexpected = 24;
			auto N = model.GetNBandsCritical();
			Assert::IsTrue(N == Nexpected);
			ArrayXf frequencies = model.GetCenterFrequencies();
			outputLog << "Center frequencies: " << frequencies << std::endl;
			Assert::IsTrue(frequencies.size() == Nexpected);
		}

		TEST_METHOD(InverseTest)
		{
			CriticalBands model;
			model.Initialize();
			auto N = model.GetNBandsCritical();
			auto c = model.GetCoefficients();
			Eigen::ArrayXf input(c.NBands), output(N);
			input.setOnes();
			model.Process(input, output);
			Assert::IsTrue((output - 1.f).abs2().sum() == 0);
			input.setZero();
			model.Inverse<float>(output, input);
			Assert::IsTrue((input - 1.f).abs2().sum() == 0);
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
			c.NChannelsLoopback = 2;
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

	TEST_CLASS(DetectVoiceActivityTest)
	{

		TEST_METHOD(Interface)
		{
			outputLog << "Running VoiceActivityDetectionTest->Interface.\n";
			DetectVoiceActivation VAD;
			auto c = VAD.GetCoefficients();
			ArrayXXcf input(c.NBands, c.NChannels);
			auto output = true;
			input.setRandom();
			auto flag = AlgorithmInterfaceTest<DetectVoiceActivation>(input, output);
			Assert::IsTrue(flag);
		}

		TEST_METHOD(ActivityDetect)
		{
			outputLog << "Running VoiceActivityDetectionTest->ActivityDetect.\n";
			DetectVoiceActivation VAD;
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