#pragma once
#include "Newton_polynomial.h"


struct Cubic_spline
{
	int n;
	double a, b;
	double(*p_func)(double x);

	double h;
	double* mas_x;	//значения x_i
	double* mas_a;	//коэффициенты сплайнов
	double* mas_b;
	double* mas_c;
	double* mas_d;

	double max_d = 0;

	Cubic_spline(int n, double a, double b, double(*p_func)(double x)) :n(n), a(a), b(b), p_func(p_func)
	{
		assert(a <= b);

		h = (b - a) / (n - 1);

		mas_x = new double[n];
		mas_a = new double[n];
		double curr_x = a;
		for (int i = 0; i < n; ++i)
		{
			mas_x[i] = curr_x;
			mas_a[i] = p_func(mas_x[i]);
			curr_x += h;
		}
		mas_b = new double[n];
		mas_c = new double[n];
		mas_d = new double[n];

		mas_c[0] = 0;		//естественные граничные условия
		mas_c[n - 1] = 0;

		for (int i = 1; i < n - 1; ++i)
		{
			mas_c[i] = (mas_a[i + 1] - 2 * mas_a[i] + mas_a[i - 1]) - (mas_c[i - 1] / 4);
		}
		for (int i = n - 2; i > 1; --i)
		{
			mas_c[i] -= mas_c[i + 1];
			mas_c[i] *= 4. / 15;
		}
		mas_c[1] -= mas_c[2];
		mas_c[1] /= 4.;
		for (int i = 0; i < n; ++i)
		{
			mas_c[i] *= 3. / (h*h);
		}

		for (int i = 1; i < n; ++i)
		{
			mas_b[i] = (mas_a[i] - mas_a[i - 1]) / h + (2 * mas_c[i] + mas_c[i - 1]) * h / 3;
			mas_d[i] = (mas_c[i] - mas_c[i - 1]) / (3 * h);
		}


		char fileName[30];
		sprintf_s(fileName, "Cube_spline_%d.txt", n);
		output(fileName);
	}

	double answer(double x)
	{
		if (x >= b)
		{
			return mas_a[n - 1];
		}
		else
		{
			int l = 0;
			while (true)
			{
				++l;
				if (mas_x[l] > x)
				{
					--l;
					break;
				}
				if (l == n - 1)
					break;
			}
			++l;
			return mas_a[l] + mas_b[l] * (x - mas_x[l]) + mas_c[l] * (x - mas_x[l])*(x - mas_x[l]) +
				mas_d[l] * (x - mas_x[l])*(x - mas_x[l])*(x - mas_x[l]);
		}
	}


	void output(char fileName[])
	{
		ofstream fout(fileName);
		fout << fixed << setprecision(10);
		double x = a, h_to_paint = (b - a) / 999, d, ans;
		int i = 0;
		while (x <= b)
		{
			ans = answer(x);
			fout << x << '\t' << ans << '\n';
			d = abs(func_values[i] - ans);
			if (d > max_d)
				max_d = d;
			x += h_to_paint;
			++i;
		}
		fout.close();
	}

};