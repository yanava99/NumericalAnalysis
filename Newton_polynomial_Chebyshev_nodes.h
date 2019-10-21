#pragma once
#include "Newton_polynomial.h"



struct Newton_Chebyshev : public Newton_polynomial
{

	Newton_Chebyshev(int n, double a, double b, double(*p_func)(double x)) : Newton_polynomial(n, a, b, p_func)
	{
		assert(a <= b);
		mas_x = new double[n];
		for (int i = 1; i <= n; ++i)
		{
			mas_x[n - i] = (a + b) / 2 + (b - a)*cos((2 * i - 1)*PI / (2 * n)) / 2;
		}

		matr = new double*[n];
		matr[0] = new double[n];
		for (int i = 0; i < n; ++i)
		{
			matr[0][i] = p_func(mas_x[i]);
		}
		for (int i = 1; i < n; ++i)
		{
			matr[i] = new double[n - i];
			for (int j = 0; j < n - i; ++j)
			{
				matr[i][j] = (matr[i - 1][j + 1] - matr[i - 1][j]) / (mas_x[j + i] - mas_x[j]);
			}
		}

		char fileName[30];
		sprintf_s(fileName, "Newton_Chebyshev_%d.txt", n);
		output(fileName);
	}

};