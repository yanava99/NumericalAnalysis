#pragma once
#include "Newton_polynomial.h"



struct Newton_equidistant : public Newton_polynomial
{

	Newton_equidistant(int n, double a, double b, double(*p_func)(double x)) : Newton_polynomial(n, a, b, p_func)
	{
		assert(a <= b);
		double h = (b - a) / (n - 1);

		mas_x = new double[n];
		mas_x[0] = a;
		for (int i = 1; i < n; ++i)
		{
			mas_x[i] = mas_x[i - 1] + h;
		}

		matr = new double*[n];
		matr[0] = new double[n];
		for (int i = 0; i < n; ++i)
		{
			matr[0][i] = p_func(mas_x[i]);
		}
		double curr_h;
		for (int i = 1; i < n; ++i)
		{
			matr[i] = new double[n - i];
			curr_h = h * i;
			for (int j = 0; j < n - i; ++j)
			{
				matr[i][j] = (matr[i - 1][j + 1] - matr[i - 1][j]) / curr_h;
			}
		}

		char fileName[30];
		sprintf_s(fileName, "Newton_equidistant_%d.txt", n);
		output(fileName);
	}

};