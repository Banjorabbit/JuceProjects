#include "unittest.h"
#include "../Utilities/pffft.h"

namespace EigenTests
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
#if NDEBUG

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
			Assert::IsTrue(error*error < 1e-8f);

		}
#endif // #if NDEBUG
	};
}
