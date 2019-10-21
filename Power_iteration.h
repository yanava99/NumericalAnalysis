#pragma once

#include "Matrix.h"
#include "My_type_matrix.h"
#include "Gradient_descent.h"
#include "Iterative_method.h"

template<class Matr>
struct Power_iter_max
{
	Matr A;
	double eigenvalue;
	Vector eigenvector;

	Power_iter_max(Matr A) : A(A)
	{
		Vector w(A.n, 1);
		do
		{
			eigenvector = w * (1.0 / w.norm());
			w = A.mult_transpose(A * eigenvector);
			eigenvalue = w.scalar(eigenvector) / eigenvector.scalar(eigenvector);
		} while ((w - eigenvector * eigenvalue).norm() >= EPS);
	}

	double norm_exp()
	{
		return (A.mult_transpose(A*eigenvector) - eigenvector * eigenvalue).norm();
	}

	void showEigenvector()
	{
		eigenvector.show();
	}

	void showA()
	{
		A.show();
	}
};


struct Power_iter_min
{
	Matrix A;
	double eigenvalue;
	Vector eigenvector;

	Power_iter_min(Matrix A) : A(A)
	{
		Vector w;
		eigenvector = Vector(A.n, 1);
		do
		{
			Gradient_descent<Matrix> gd1(A.transpose(), eigenvector);
			Gradient_descent<Matrix> gd2(A, gd1.x);
			w = gd2.x;
			eigenvalue = w.scalar(eigenvector) / eigenvector.scalar(eigenvector);
			eigenvector = w * (1.0 / w.norm());
		} while ((w - eigenvector * eigenvalue).norm() >= EPS_POWER);
		eigenvalue = 1.0 / eigenvalue;
	}

	double norm_exp()
	{
		return (A.mult_transpose(A * eigenvector) - eigenvector * eigenvalue).norm();
	}

	void showEigenvector()
	{
		eigenvector.show();
	}

	void showA()
	{
		A.show();
	}
};