#pragma once
#include "Matrix.h"

const double EPS_Gauss = 0.000000000000001;

struct Gaussian_elimination
{
	Matrix A;
	Vector b;


	Gaussian_elimination(Matrix A, Vector b) : A(A), b(b)
	{
		if (A.n != b.n)
			throw invalid_argument("Dimensions aren't compatible.");
		else
		{
			findX();
		}
	}

	void findX()
	{
		//прямой ход метода Гаусса
		int index_max = 0;
		for (int i = 0; i < A.n; ++i)
		{
			//выбор главного элемента
			index_max = i;
			for (int j = i; j < A.n; ++j)
				if (abs(A.A[j][i]) > abs(A.A[index_max][i]))
					index_max = j;

			if (abs(A.A[index_max][i]) < EPS_Gauss)
				throw invalid_argument("Matrix is invalid.");

			swap(A.A[i], A.A[index_max]);
			swap(b.x[i], b.x[index_max]);

			for (int j = i + 1; j < A.n; ++j)
				A.A[i][j] /= A.A[i][i];
			b.x[i] /= A.A[i][i];
			for (int j = i + 1; j < A.n; ++j)
			{
				for (int k = i + 1; k < A.n; ++k)
					A.A[j][k] -= A.A[i][k] * A.A[j][i];
				b.x[j] -= b.x[i] * A.A[j][i];
			}
		}

		//обратный ход метода Гаусса
		for (int j = A.n - 1; j > -1; --j)
			for (int i = j - 1; i > -1; --i)
				b.x[i] -= b.x[j] * A.A[i][j];

	}

	void showX()
	{
		b.show();
	}
};