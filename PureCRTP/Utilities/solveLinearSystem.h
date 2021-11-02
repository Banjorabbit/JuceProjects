#pragma once

// solve linear system of equations Ax = b and save result in b
void solveLinearSystemLDLT(float *A, float *b, int N, int Nmic);
void solveLinearSystemLLT(float *A, float *b, int N, int Nmic);
void solveLinearSystemQR(float *A, float *b, int N, int Nmic);
void solveLinearSystemQRCol(float *A, float *b, int N, int Nmic);
void solveLinearSystemQRFull(float *A, float *b, int N, int Nmic);
void solveLinearSystemLU(float *A, float *b, int N, int Nmic);
void solveLinearSystemLUFull(float *A, float *b, int N, int Nmic);
void solveLinearSystemSVD(float *A, float *b, int N, int Nmic);

