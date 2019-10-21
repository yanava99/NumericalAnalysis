#pragma once
#include <math.h>
#include <assert.h>



struct Trapezoidal_rule
{
	double(*f)(double x);
	double(*p)(double x);
	double a, b;

	double eps;
	int n;
	double ans_prev;
	double ans;


	Trapezoidal_rule(double(*f)(double x), double(*p)(double x), double a, double b, double eps)
		:f(f), p(p), a(a), b(b), eps(eps), n(2)
	{
		double h = (b - a) / n;
		ans = (p(a)*f(a) + p(b)*f(b)) / 2;
		double x = a;
		for (int i = 1; i < n; ++i)
		{
			x += h;
			ans += p(x)*f(x);
		}
		ans *= h;
		do
		{
			ans_prev = ans;
			n *= 2;
			h = (b - a) / n;
			ans = (p(a)*f(a) + p(b)*f(b)) / 2;
			x = a;
			for (int i = 1; i < n; ++i)
			{
				x += h;
				ans += p(x)*f(x);
			}
			ans *= h;
		} while (abs(ans - ans_prev) >= 3 * eps);
	}
};