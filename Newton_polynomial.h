#pragma once
#include <assert.h>
#include <math.h>
#include <fstream>
#include <iomanip>

using namespace std;



const double PI = 3.141592653589793238463;

double* func_values = new double[1000];




struct Newton_polynomial
{
	int n;
	double a, b;
	double(*p_func)(double x);

	double* mas_x;		//значение x_i
	double** matr;		//разделённые разности

	double max_d = 0;

	Newton_polynomial(int n, double a, double b, double(*p_func)(double x)) :
		n(n), a(a), b(b), p_func(p_func) {}


	//взятие первых разделённых разностей

	//double answer(double x)
	//{
	//	double ans = matr[0][0], curr_ans = 1;
	//	for (int i = 1; i < n; ++i)
	//	{
	//		curr_ans *= x - mas_x[i - 1];
	//		ans += curr_ans * matr[i][0];
	//	}
	//	return ans;
	//}


	//взятие разделённых разностей, близких к точке x

	double answer(double x)
	{
		if (x >= b)
		{
			return matr[0][n - 1];
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
			double curr_ans = 1;
			double ans = matr[0][l];
			for (int i = 0; i < n - l - 1; ++i)
			{
				curr_ans *= x - mas_x[l + i];
				ans += curr_ans * matr[i + 1][l];
			}
			for (int i = n - l; i < n; ++i)
			{
				curr_ans *= x - mas_x[n - i - 1];
				ans += curr_ans * matr[i][n - i - 1];
			}
			return ans;
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