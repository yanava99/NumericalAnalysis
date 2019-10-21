#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

const double EPS = 0.000001;
const double EPS_POWER = 0.00000000001;
unsigned int mult_count = 0;
unsigned int effective_mult_count = 0;

struct Vector
{
	int n;
	vector<double> x;

	Vector() : n(0) {}

	Vector(int n) : n(n)
	{
		for (int i = 0; i < n; ++i)
			x.push_back(0);
	}

	Vector(int n, double number) : n(n)
	{
		for (int i = 0; i < n; ++i)
			x.push_back(number);
	}

	Vector(const Vector& toCopy)
	{
		n = toCopy.n;
		x = toCopy.x;
	}

	Vector(vector<double> v)
	{
		n = v.size();
		x = v;
	}

	Vector(double x_0, double x_1, double x_2, double x_3) : n(4)
	{
		x.push_back(x_0);
		x.push_back(x_1);
		x.push_back(x_2);
		x.push_back(x_3);
	}

	Vector operator+(const Vector& other) const
	{
		Vector res(n);
		if (n != other.n)
			throw invalid_argument("Dimensions aren't compatible.");
		else
		{
			for (int i = 0; i < n; ++i)
				res.x[i] = x[i] + other.x[i];
		}
		return res;
	}

	Vector operator-(const Vector& other) const
	{
		Vector res(n);
		if (n != other.n)
			throw invalid_argument("Dimensions aren't compatible.");
		else
		{
			for (int i = 0; i < n; ++i)
				res.x[i] = x[i] - other.x[i];
		}
		return res;
	}

	Vector operator-() const
	{
		Vector res(n);
		for (int i = 0; i < n; ++i)
			res.x[i] = -x[i];
		return res;
	}

	double scalar(const Vector& other) const
	{
		double res = 0;
		if (n != other.n)
			throw invalid_argument("Dimensions aren't compatible.");
		else
		{
			for (int i = 0; i < n; ++i)
				res += x[i] * other.x[i];
		}
		return res;
	}

	Vector operator*(double number) const
	{
		Vector res(n);
		for (int i = 0; i < n; ++i)
			res.x[i] = x[i] * number;
		return res;
	}

	double norm() const
	{
		double max = fabs(x[0]);
		for (int i = 1; i < n; ++i)
			if (fabs(x[i]) > max)
				max = fabs(x[i]);

		return max;
	}

	void show()
	{
		for (int i = 0; i < n; ++i)
			cout << x[i] << '\n';
	}
};
