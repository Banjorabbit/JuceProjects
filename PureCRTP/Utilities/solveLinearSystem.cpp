#include "solveLinearSystem.h"

#define EIGEN_NO_DEBUG 
#include "../eigen/Eigen/Dense"
using namespace Eigen;

void solveLinearSystemLDLT(float *A, float *b, int N, int Nmic)
{
	Map<MatrixXcf>((std::complex<float>*)(b), N, Nmic) = Map<MatrixXcf>((std::complex<float>*)(A), N, N).ldlt().solve(Map<MatrixXcf>((std::complex<float>*)(b), N, Nmic));
}

void solveLinearSystemLLT(float *A, float *b, int N, int Nmic)
{
	Map<MatrixXcf>((std::complex<float>*)(b), N, Nmic) = Map<MatrixXcf>((std::complex<float>*)(A), N, N).llt().solve(Map<MatrixXcf>((std::complex<float>*)(b), N, Nmic));
}

void solveLinearSystemQR(float *A, float *b, int N, int Nmic)
{
	Map<MatrixXcf>((std::complex<float>*)(b), N, Nmic) = Map<MatrixXcf>((std::complex<float>*)(A), N, N).householderQr().solve(Map<MatrixXcf>((std::complex<float>*)(b), N, Nmic));
}

void solveLinearSystemQRCol(float *A, float *b, int N, int Nmic)
{
	Map<MatrixXcf>((std::complex<float>*)(b), N, Nmic) = Map<MatrixXcf>((std::complex<float>*)(A), N, N).colPivHouseholderQr().solve(Map<MatrixXcf>((std::complex<float>*)(b), N, Nmic));
}

void solveLinearSystemQRFull(float *A, float *b, int N, int Nmic)
{
	Map<MatrixXcf>((std::complex<float>*)(b), N, Nmic) = Map<MatrixXcf>((std::complex<float>*)(A), N, N).fullPivHouseholderQr().solve(Map<MatrixXcf>((std::complex<float>*)(b), N, Nmic));
}

void solveLinearSystemLU(float *A, float *b, int N, int Nmic)
{
	Map<MatrixXcf>((std::complex<float>*)(b), N, Nmic) = Map<MatrixXcf>((std::complex<float>*)(A), N, N).partialPivLu().solve(Map<MatrixXcf>((std::complex<float>*)(b), N, Nmic));
}

void solveLinearSystemLUFull(float *A, float *b, int N, int Nmic)
{
	Map<MatrixXcf>((std::complex<float>*)(b), N, Nmic) = Map<MatrixXcf>((std::complex<float>*)(A), N, N).fullPivLu().solve(Map<MatrixXcf>((std::complex<float>*)(b), N, Nmic));
}

void solveLinearSystemSVD(float *A, float *b, int N, int Nmic)
{
	Map<MatrixXcf>((std::complex<float>*)(b), N, Nmic) = Map<MatrixXcf>((std::complex<float>*)(A), N, N).bdcSvd(ComputeThinU | ComputeThinV).solve(Map<MatrixXcf>((std::complex<float>*)(b), N, Nmic));
}
