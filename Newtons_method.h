#pragma once
#include "Gaussian_elimination.h"
#include <fstream>
#include <iomanip>
#include <cmath>


const double EPS_NEWTON = 0.0000000001;

struct Newtons_method_for_my_function
{
	Vector x;

	Newtons_method_for_my_function(Vector x)
	{
		if (x.n != 4)
			throw invalid_argument("Dimensions aren't compatible.");
		this->x = x;

		findX();
	}

	Newtons_method_for_my_function(double A0, double A1, double x0, double x1)
	{
		x = Vector(4);
		x.x[0] = A0;
		x.x[1] = A1;
		x.x[2] = x0;
		x.x[3] = x1;

		findX();
	}


	Matrix J(Vector v)
	{
		if (x.n != 4)
			throw invalid_argument("Dimensions aren't compatible.");

		Matrix J = Matrix(4);

		J.A[0][0] = 1.;
		J.A[0][1] = 1.;
		J.A[0][2] = 0.;
		J.A[0][3] = 0.;

		J.A[1][0] = v.x[2];
		J.A[1][1] = v.x[3];
		J.A[1][2] = v.x[0];
		J.A[1][3] = v.x[1];

		J.A[2][0] = v.x[2] * v.x[2];
		J.A[2][1] = v.x[3] * v.x[3];
		J.A[2][2] = 2. * v.x[0] * v.x[2];
		J.A[2][3] = 2. * v.x[1] * v.x[3];

		J.A[3][0] = v.x[2] * v.x[2] * v.x[2];
		J.A[3][1] = v.x[3] * v.x[3] * v.x[3];
		J.A[3][2] = 3. * v.x[0] * v.x[2] * v.x[2];
		J.A[3][3] = 3. * v.x[1] * v.x[3] * v.x[3];

		return J;
	}

	void findX()
	{
		ofstream fout("output.txt");
		fout << setprecision(16) << fixed;
		setlocale(LC_ALL, "Russian");
		int k = 1;
		double norm = F(x).norm();
		do
		{
			Vector f = -F(x);
			Gaussian_elimination gauss(J(x), f);
			x = x + gauss.b;
			norm = F(x).norm();
			fout /*<< k */ << ' ' << norm << '\n';
			++k;
		} while (norm > EPS_NEWTON);
	}

	Vector F(Vector x)
	{
		if (x.n != 4)
			throw invalid_argument("Dimensions aren't compatible.");

		Vector f(
			x.x[0] + x.x[1] - 4. / 3,
			x.x[0] * x.x[2] + x.x[1] * x.x[3],
			x.x[0] * x.x[2] * x.x[2] + x.x[1] * x.x[3] * x.x[3] - 4. / 15,
			x.x[0] * x.x[2] * x.x[2] * x.x[2] + x.x[1] * x.x[3] * x.x[3] * x.x[3]);

		return f;
	}
};