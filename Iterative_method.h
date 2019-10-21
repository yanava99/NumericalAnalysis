#pragma once
#include "Matrix.h"
#include "My_type_matrix.h"


template<class Matr>
struct Iterative_method
{
	Matr A;
	Vector b, x;

	Iterative_method(string fileName)
	{
		A.fromFile(fileName);
		setB();
		findX();
	}

	Iterative_method(Matr A) : A(A)
	{
		setB();
		findX();
	}

	Iterative_method(Matr A, bool generate_data) : A(A)
	{
		setB();
		if (generate_data)
			findX_data();
		else
			findX();
	}

	Iterative_method(Matr A, Vector b) : A(A), b(b)
	{
		if (A.n != b.n)
			throw invalid_argument("Dimensions aren't compatible.");
		else
		{
			findX();
		}
	}

	void setB()
	{
		Vector e(A.n, 1);
		b = A * e;
	}

	void findX()
	{
		x = Vector(A.n);
		Vector temp = A * x - b;
		do
		{
			x = x - (A.mult_transpose(temp)*(1.0 / A.norm_of_transpose_mult_A));
			temp = A * x - b;
		} while (temp.norm() >= EPS);
	}

	void findX_data()
	{
		ofstream out("iter_method_data.dat");
		x = Vector(A.n);
		Vector temp = A * x - b;
		int i = 0;
		double norm;
		do
		{
			x = x - (A.mult_transpose(temp)*(1.0 / A.norm_of_transpose_mult_A));
			temp = A * x - b;
			norm = temp.norm();
			++i;
			if (i % 20000 == 0)
			{
				out << i << ' ' << norm << '\n';
			}
		} while (norm >= EPS);
		out.close();
	}

	void showX()
	{
		x.show();
	}

	void showA()
	{
		A.show();
	}

	void showB()
	{
		b.show();
	}
};