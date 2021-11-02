#include "unittest.h"
#include "../Utilities/pffft.h"
#include "../Utilities/HelperFunctions.h"

namespace EigenTests
{
	TEST_CLASS(EigenMethodsTest)
	{
		TEST_METHOD(HeadZeroTest)
		{
			ArrayXf x(4);
			x << 1.f, 2.f, 3.f, 4.f;
			int n = 0;
			outputLog << "x.head(0).sum(): " << x.head(n).sum() << std::endl;
		}

		TEST_METHOD(Index1DimensionIn2DimensionalArrayTest)
		{
			ArrayXXf A(3, 2);
			A << 1.f, 2.f, 3.f, 4.f, 5.f, 6.f;
			outputLog << "A: " << A << std::endl;
			outputLog << "A(0): " << A(0) << ", A(1): " << A(1) << ", A(2): " << A(2) << ", A(3): " << A(3) << ", A(4): " << A(4) << ", A(5): " << A(5) << std::endl;
			outputLog << "A(0,0): " << A(0,0) << ", A(1,0): " << A(1,0) << ", A(2,0): " << A(2,0) << ", A(0,1): " << A(0,1) << ", A(1,1): " << A(1,1) << ", A(2,1): " << A(2,1) << std::endl;
		}

		TEST_METHOD(MultiplyWithBool)
		{
			ArrayXf A(4);
			A << 1.f, 2.f, 3.f, 4.f;
			Array<bool, Dynamic, 1> B(4);
			B << true, false, true, false;
			outputLog << A * B.cast<float>() << std::endl;
		}

		TEST_METHOD(ZeroBeforeResize)
		{
			ArrayXXf A;
			A.setZero(); // set to zero even though it does not have a size.
		}

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

		TEST_METHOD(MapTest)
		{
			int size = 10;
			ArrayXf A(size), B(size);
			float ref[] = { 0,1,2,3,4,5,6,7,8,9 };
			A.setOnes();
			B.setZero();
			Map<ArrayXf>(A.data(),size) = Map<ArrayXf>(B.data(), size);
			Map<ArrayXf>(B.data(),size) = Map<ArrayXf>(ref, size);
			

			ref[0] = 10;
			B(1) = 11;
			A(2) = 12;
			outputLog << "A: " << A << "\n";
			outputLog << "B: " << B << std::endl;
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
		TEST_METHOD(AngleTest)
		{
			ArrayXf test(3);
			test.setZero();
			outputLog << "Angle: " << test.arg() << std::endl;
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
			ArrayXcf output(N / 2 + 1);

			pffft_transform_ordered(setup, &input(0), &input(0), NULL, directionFFT);
			output(0) = input(0);
			output.tail(1) = input(1);
			std::memcpy(&output.real()(1), &input(2), (N - 2) * sizeof(float));
			outputLog << output << std::endl;

		}

	};

	TEST_CLASS(Efficiency)
	{
		// calculate the diagonal elements of the inverse matrix. This can be done in two ways:
		// 1) Calculate the full inverse matrix and then take the diagonal elements
		// 2) Calculate the determinant of a subset of the full matrix and divide by the determinant of the full matrix
		//
		// This test shows the two methods and compares the execution speed
		TEST_METHOD(InverseDiagonalTest)
		{
			double durationLoop = 1e10;
			double durationInverse = 1e10;
			int NChannels = 4; // execution times depend heavily on this value
			MatrixXcf EigenVectors(NChannels, NChannels);
			EigenVectors.setRandom();
			VectorXcf invValues(NChannels);
			MatrixXcf newMatrix(NChannels - 1, NChannels - 1);
			VectorXcf EigenInv(NChannels);

			auto start = std::chrono::steady_clock::now();
			for (auto i = 0; i < NChannels; i++)
			{
				int ai = 0;
				for (auto a = 0; a < NChannels; a++)
				{
					if (a != i)
					{
						int bi = 0;
						for (auto b = 0; b < NChannels; b++)
						{
							if (b != i)
							{
								newMatrix(ai, bi) = EigenVectors(a, b);
								bi++;
							}
						}
						ai++;
					}
				}
				//outputLog << "newMatrix: " << newMatrix << std::endl; // uncomment this to print newMatrix in each for-loop
				invValues(i) = newMatrix.determinant();
			}
			invValues /= EigenVectors.determinant();
			auto end = std::chrono::steady_clock::now();
			auto time = std::chrono::duration<double, std::micro>(end - start).count();
			durationLoop = std::min(durationLoop, time);
			
			start = std::chrono::steady_clock::now();
			EigenInv = EigenVectors.inverse();
			end = std::chrono::steady_clock::now();
			time = std::chrono::duration<double, std::micro>(end - start).count();
			durationInverse = std::min(durationInverse, time);

			// only calculate one value for reference
			double durationSingle = 1e10;
			start = std::chrono::steady_clock::now();
			std::complex<float> EigenSingleInv = EigenVectors.block(1, 0, NChannels - 1, NChannels - 1).determinant() / EigenVectors.determinant();
			end = std::chrono::steady_clock::now();
			time = std::chrono::duration<double, std::micro>(end - start).count();
			durationSingle = std::min(durationSingle, time);

			// Fixed size matrix
			Matrix4cf EigenVectorsFix;
			EigenVectorsFix.setRandom();
			Matrix4cf EigenInvFix;
			
			double durationFix = 1e10;
			start = std::chrono::steady_clock::now();
			std::complex<float> EigenFixInv = EigenVectorsFix.block(1, 0, NChannels - 1, NChannels - 1).determinant() / EigenVectorsFix.determinant();
			end = std::chrono::steady_clock::now();
			time = std::chrono::duration<double, std::micro>(end - start).count();
			durationFix = std::min(durationSingle, time);

			double durationInvFix = 1e10;
			start = std::chrono::steady_clock::now();
			EigenInvFix = EigenVectorsFix.inverse();
			end = std::chrono::steady_clock::now();
			time = std::chrono::duration<double, std::micro>(end - start).count();
			durationInvFix = std::min(durationSingle, time);

			outputLog << "Values: " << invValues << std::endl;
			outputLog << "Inverse: " << EigenInv << std::endl;

			outputLog << "Duration loop: " << durationLoop << std::endl;
			outputLog << "Duration inverse: " << durationInverse << std::endl;
			outputLog << "Duration single: " << durationSingle << std::endl;

			outputLog << "Duration fixed single: " << durationFix << std::endl;
			outputLog << "Duration fixed inverse: " << durationInvFix << std::endl;
		}

//#if NDEBUG

		TEST_METHOD(AssignMult)
		{
			int nChannels = 257;
			int size = 257;
			int index = 0;
			int length = 257;
			Eigen::ArrayXXf Output(nChannels, nChannels);
			Eigen::ArrayXXf Input(nChannels, nChannels);
			Eigen::ArrayXXf Input2(nChannels, nChannels);
			Output.setZero();
			Input.setRandom();
			Input2.setRandom();

			auto start = std::chrono::steady_clock::now();
			//for (auto channel = 0;channel < size;channel++)
			//{
			//	Output.block(index,channel, length, 1) = Input.block(index, channel, length, 1) * Input2.block(index, channel, length, 1);
			//	//Output.col(channel) = Input.col(channel) * Input2.col(channel);
			//}
			Output.block(index, 0, length, size) =Input.block(index, 0, length, size) * Input2.block(index, 0, length, size);
			//Output = Input * Input2;
			auto end = std::chrono::steady_clock::now();
			auto time1 = static_cast<double>(std::chrono::duration<double, std::micro>(end - start).count());
			decltype(Output) res1 = Output;

			start = std::chrono::steady_clock::now();
			for (auto channel = 0;channel < size;channel++)
			{
				float *p1 = &Output(index, channel);
				float *p2 = &Input(index, channel);
				float *p3 = &Input2(index, channel);
				Map<ArrayXf>(p1, length) = Map<ArrayXf>(p2, length) * Map<ArrayXf>(p3, length);
			}
			end = std::chrono::steady_clock::now();
			auto time2 = static_cast<double>(std::chrono::duration<double, std::micro>(end - start).count());
			decltype(Output) res2 = Output;

			start = std::chrono::steady_clock::now();
			for (auto channel = 0;channel < size;channel++)
			{
				float *p1 = &Output(index, channel);
				float *p2 = &Input(index, channel);
				float *p3 = &Input2(index, channel);
				assignMultEigen(p1, p2, p3, length);
			}
			//float *p1 = &Output(index, 0);
			//float *p2 = &Input(index, 0);
			//float *p3 = &Input2(index, 0);
			//assignMultEigen(p1, p2, p3, length*size);
			end = std::chrono::steady_clock::now();
			auto time3 = static_cast<double>(std::chrono::duration<double, std::micro>(end - start).count());
			decltype(Output) res3 = Output;


			outputLog << "Time 1: " << time1 << "\n";
			outputLog << "Time 2: " << time2 << "\n";
			outputLog << "Time 3: " << time3 << "\n";
			outputLog << "Error2: " << (res2-res1).abs2().sum() << "\n";
			outputLog << "Error3: " << (res3-res1).abs2().sum() << std::endl;
		}

		TEST_METHOD(AssignForLoopVsEigen)
		{
			// put circular loopback buffer into non-circular affine matrix
			int filterLength = 120;
			int orderAffine = 5;
			int bufferLength = filterLength + orderAffine - 1;
			int circCounter = 32;
			int affineCounter = 2;
			int nBands = 257;

			Eigen::ArrayXXcf buffersLoopback(bufferLength, nBands);
			buffersLoopback.setRandom();
			std::vector<Eigen::ArrayXXcf> bufferAffine;
			bufferAffine.resize(nBands);
			for (auto &buffer : bufferAffine) { buffer.resize(filterLength, orderAffine); }

			// output test data
			std::vector<Eigen::ArrayXXcf> outTest;
			outTest.resize(nBands);
			for (auto &buffer : outTest) { buffer.resize(filterLength, orderAffine); }
			
			const int bufferLength1 = std::min(bufferLength - circCounter, filterLength);
			const int bufferLength2 = filterLength - bufferLength1;

			auto start = std::chrono::steady_clock::now();
			for (auto ibin = 0; ibin < nBands; ibin++)
			{
				bufferAffine[ibin].col(affineCounter).head(bufferLength1) = buffersLoopback.col(ibin).segment(circCounter, bufferLength1);
				bufferAffine[ibin].col(affineCounter).tail(bufferLength2) = buffersLoopback.col(ibin).head(bufferLength2);
			}
			auto end = std::chrono::steady_clock::now();
			auto duration1 = static_cast<double>(std::chrono::duration<double, std::micro>(end - start).count());

			// save output
			for (auto ibin = 0; ibin < nBands; ibin++)
			{
				outTest[ibin] = bufferAffine[ibin];
				bufferAffine[ibin] = 0;
			}

			start = std::chrono::steady_clock::now();
			std::complex<float> *bp2 = buffersLoopback.data();
			for (auto ibin = 0; ibin < nBands; ibin++, bp2 += bufferLength)
			{
				std::complex<float> *bp1 = bufferAffine[ibin].col(affineCounter).data();
				Eigen::Map<ArrayXcf>(bp1, bufferLength1) = Eigen::Map<ArrayXcf>(bp2 + circCounter, bufferLength1);
				Eigen::Map<ArrayXcf>(bp1 + bufferLength1, bufferLength2) = Eigen::Map<ArrayXcf>(bp2, bufferLength2);
			}
			end = std::chrono::steady_clock::now();
			auto duration2 = static_cast<double>(std::chrono::duration<double, std::micro>(end - start).count());

			float error2 = 0;
			for (auto ibin = 0; ibin < nBands; ibin++)
			{
				error2 += (bufferAffine[ibin].col(affineCounter) - outTest[ibin].col(affineCounter)).abs2().mean();
				bufferAffine[ibin] = 0;
			}
			outputLog << "Error2: " << error2 << "\n";

			start = std::chrono::steady_clock::now();
			
			for (auto ibin = 0; ibin < nBands; ibin++)
			{
				std::complex<float> *bp1 = bufferAffine[ibin].col(affineCounter).data();
				bp2 = buffersLoopback.col(ibin).data();
				assignCircularEigen(bp1, bp2, circCounter, bufferLength1, bufferLength2);
			}
			end = std::chrono::steady_clock::now();
			auto duration3 = static_cast<double>(std::chrono::duration<double, std::micro>(end - start).count());

			float error3 = 0;
			for (auto ibin = 0; ibin < nBands; ibin++)
			{
				error3 += (bufferAffine[ibin].col(affineCounter) - outTest[ibin].col(affineCounter)).abs2().mean();
			}
			outputLog << "Error3: " << error3 << "\n";

			outputLog << "Time 1: " << duration1 << "\n";
			outputLog << "Time 2: " << duration2 << "\n";
			outputLog << "Time 3: " << duration3 << "\n";
			outputLog << std::endl;
		}
		TEST_METHOD(EigenValues)
		{
			srand((unsigned int)time(0));
			EigenSolver<MatrixXf> es;
			int N = 20;
			es = EigenSolver<MatrixXf>::EigenSolver(N);

			MatrixXf mat(N, N);
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
			Assert::IsTrue(error*error < 1e-8f);

		}
//#endif // #if NDEBUG
	};
}
