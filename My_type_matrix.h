#pragma once
#include "Matrix.h"


//this type was specified in initial condition

//   0   0   0  10   0   1
//   0   0  10   0   2 -10
//   0  10   0   3 -10   0
//  10   0   4 -10   0   0
//   0   5 -10   0   0   0
//   6 -10   0   0   0   0



struct My_type_matrix : public Matrix
{
	My_type_matrix(int n)
	{
		if (n < 4)
			throw invalid_argument("Can't create matrix with negative dimension.");
		else
		{
			this->n = n;
			vector<double> temp;
			for (int i = 0; i < n; ++i)
			{
				for (int j = 0; j < n; ++j)
				{
					if (j == n - i - 1)
						temp.push_back(i + 1);
					else if (j == n - i - 3)
						temp.push_back(10);
					else if (j == n - i)
						temp.push_back(-10);
					else
						temp.push_back(0);
				}
				A.push_back(temp);
				temp.clear();
			}
			set_norm_of_transpose_mult_A();
		}
	}

	void fromFile(string fileName)
	{
		throw invalid_argument("This method is not available.");
	}

	My_type_matrix(const My_type_matrix& other)
	{
		this->n = other.n;
		this->A = other.A;
		this->norm_of_transpose_mult_A = other.norm_of_transpose_mult_A;
	}

	Vector operator*(const Vector& x) const
	{
		Vector res(n);
		if (n != x.n)
			throw invalid_argument("Dimensions aren't compatible.");
		else
		{
			res.x[0] = 10 * x.x[n - 3] + x.x[n - 1];
			for (int i = 1; i < n - 2; ++i)
				res.x[i] = 10 * x.x[n - i - 3] + (i + 1) * x.x[n - i - 1] - 10 * x.x[n - i];
			res.x[n - 2] = (n - 1) * x.x[1] - 10 * x.x[2];
			res.x[n - 1] = n * x.x[0] - 10 * x.x[1];
			++effective_mult_count;
		}
		return res;

	}

	double norm() const
	{
		return (n + 18);
	}

	Vector mult_transpose(const Vector& x) const
	{
		Vector res(n);
		if (n != x.n)
			throw invalid_argument("Dimensions aren't compatible.");
		else
		{
			res.x[0] = 10 * x.x[n - 3] + n * x.x[n - 1];
			for (int i = 1; i < n - 2; ++i)
				res.x[i] = 10 * x.x[n - i - 3] + (n - i) * x.x[n - i - 1] - 10 * x.x[n - i];
			res.x[n - 2] = 2 * x.x[1] - 10 * x.x[2];
			res.x[n - 1] = x.x[0] - 10 * x.x[1];
			++effective_mult_count;
		}
		return res;
	}

	void set_norm_of_transpose_mult_A()
	{
		norm_of_transpose_mult_A = n * n + 100 * n + 180;
	}
};