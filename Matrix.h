#pragma once
#include "Vector.h"
#include <fstream>
#include <string>

struct Matrix
{
	int n;
	vector<vector<double>> A;

	double norm_of_transpose_mult_A;

	Matrix() : n(0) {}

	Matrix(const Matrix& other)
	{
		this->n = other.n;
		this->A = other.A;
		this->norm_of_transpose_mult_A = other.norm_of_transpose_mult_A;
	}

	Matrix(int n) : n(n)
	{
		vector<double> temp;
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
			{
				temp.push_back(0);
			}
			A.push_back(temp);
			temp.clear();
		}
		set_norm_of_transpose_mult_A();
	}

	void fromFile(string& fileName)
	{
		ifstream in(fileName);
		in >> n;
		double inp;
		vector<double> temp;
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
			{
				in >> inp;
				temp.push_back(inp);
			}
			A.push_back(temp);
			temp.clear();
		}
		in.close();
		set_norm_of_transpose_mult_A();
	}

	Vector operator*(const Vector& x) const
	{
		Vector res(n);
		if (n != x.n)
			throw invalid_argument("Dimensions aren't compatible.");
		else
		{
			for (int i = 0; i < n; ++i)
			{
				for (int j = 0; j < n; ++j)
					res.x[i] += A[i][j] * x.x[j];
			}
		}
		++mult_count;
		return res;
	}

	double norm() const
	{
		double norm = 0, current = 0;
		for (int i = 0; i < n; ++i)
		{
			current = 0;
			for (int j = 0; j < n; ++j)
				current += abs(A[i][j]);
			if (current > norm)
				norm = current;
		}
		return norm;
	}

	void show()
	{
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
				cout << A[i][j] << ' ';
			cout << '\n';
		}
	}

	Vector mult_transpose(const Vector& x) const
	{
		Vector res(n);
		if (n != x.n)
			throw invalid_argument("Dimensions aren't compatible.");
		else
		{
			for (int i = 0; i < n; ++i)
			{
				for (int j = 0; j < n; ++j)
					res.x[i] += A[j][i] * x.x[j];
			}
			++mult_count;
		}
		return res;
	}

	void set_norm_of_transpose_mult_A()
	{
		double res = 0, cur, accum;
		Vector temp(n);
		for (int i = 0; i < n; ++i)
		{
			cur = 0;
			for (int j = 0; j < n; ++j)
			{
				accum = 0;
				for (int k = 0; k < n; ++k)
					accum += A[k][i] * A[k][j];
				cur += abs(accum);
			}
			if (cur > res)
				res = cur;
		}
		norm_of_transpose_mult_A = res;
	}

	Matrix transpose() const
	{
		Matrix res;
		res.n = n;
		vector<double> temp;
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
			{
				temp.push_back(A[j][i]);
			}
			res.A.push_back(temp);
			temp.clear();
		}
		res.set_norm_of_transpose_mult_A();
		return res;
	}
};