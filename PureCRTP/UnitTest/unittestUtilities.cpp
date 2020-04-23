#include "unittest.h"
#include "../Utilities/AudioFile.h"
#include "../Utilities/FormatHandling.h"
#include "../Utilities/pffft.h"
#include "../Utilities/ApproximationMath.h"

namespace UtilitiesTests
{

	LoggerOutput outputLog;

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