#pragma once
#include "Matrix.h"
#include "My_type_matrix.h"
#include "Vector.h"


template<class Matr>
struct Gradient_descent
{
	Matr A;
	Vector b, x;


	Gradient_descent(string fileName)
	{
		A.fromFile(fileName);
		setB();
		findX();
	}

	Gradient_descent(Matr A, bool generate_data) : A(A)
	{
		setB();
		if (generate_data)
			findX_data();
		else
			findX();
	}

	Gradient_descent(Matr A) : A(A)
	{
		setB();
		findX();
	}

	Gradient_descent(Matr A, Vector b) : A(A), b(b)
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
		x = Vector(A.n);
		Vector temp = A * x - b;
		do
		{
			Vector r = A.mult_transpose(temp);
			x = x - r * (r.scalar(r) / (A.mult_transpose(A * r)).scalar(r));
			temp = A * x - b;
		} while (temp.norm() >= EPS);
	}

	void findX_data()
	{
		ofstream out("grand_descent_data.dat");
		x = Vector(A.n);
		Vector temp = A * x - b;
		int i = 0;
		double norm;
		do
		{
			Vector r = A.mult_transpose(temp);
			x = x - r * (r.scalar(r) / (A.mult_transpose(A * r)).scalar(r));
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

	void setB()
	{
		Vector e(A.n, 1);
		b = (A * e);
	}
};