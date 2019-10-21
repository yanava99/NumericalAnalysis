#define _USE_MATH_DEFINES

#include "Gradient_descent.h"
#include "Iterative_method.h"
#include "Power_iteration.h"
#include "Newtons_method.h"
#include "Newton_polynomial_equidistant_nodes.h"
#include "Newton_polynomial_Chebyshev_nodes.h"
#include "Cubic_spline.h"
#include "Trapezoidal_rule.h"

#include <iomanip>
#include <ctime>


using namespace std;


double func(double x)
{
	return x * x + 2 * sin(10 * x);
}

double p(double x)
{
	assert(x >= -1 && x <= 1);
	return sqrt(1 - x * x);
}

int numberF = 0;

double f(double x)
{
	++numberF;
	//return pow(M_E, -x * x)*sin(x);
	return abs(pow(M_E, -x * x)*sin(x));
}


struct Program
{
	Program()
	{
		//experiment_iterative_and_gradient(10);
		//experiment_iterative_and_gradient(100);
		//experiment_iterative_and_gradient(500);

		//diagram_data(500);

		//power(10);

		//newtons();

		//interpolation();

		trapezoidal_rule();
	}

	void experiment_iterative_and_gradient(int n)
	{
		My_type_matrix A_effective(n);
		Matrix A(A_effective);

		mult_count = 0;
		unsigned int start_time = clock();
		Gradient_descent<Matrix> gd(A);
		unsigned int time_needed_gd = clock() - start_time,
			mult_count_gd = mult_count;
		//gd.showA();
		//gd.showX();

		mult_count = 0;
		start_time = clock();
		Iterative_method<Matrix> im(A);
		unsigned int time_needed_im = clock() - start_time,
			mult_count_im = mult_count;
		//im.showA();
		//im.showX();


		effective_mult_count = 0;
		start_time = clock();
		Gradient_descent<My_type_matrix> gde(A_effective);
		unsigned int time_needed_gd_effective = clock() - start_time,
			mult_count_gd_effective = effective_mult_count;
		//gde.showA();
		//gde.showX();

		effective_mult_count = 0;
		start_time = clock();
		Iterative_method<My_type_matrix> ime(A_effective);
		unsigned int time_needed_im_effective = clock() - start_time,
			mult_count_im_effective = effective_mult_count;
		//ime.showA();
		//ime.showX();

		cout << "\n\n\n---------------- Speed comparison for n = " << n << "----------------------------\n\n\n";
		cout << "Time:\n\n";
		cout << "Gradient descent: " << (double)time_needed_gd / CLOCKS_PER_SEC << '\n';
		cout << "Iterative method: " << (double)time_needed_im / CLOCKS_PER_SEC << '\n';
		cout << "Gradient descent (multiplication is implemented efficiently for a matrix of my type): " << (double)time_needed_gd_effective / CLOCKS_PER_SEC << '\n';
		cout << "Iterative method (multiplication is implemented efficiently for a matrix of my type): " << (double)time_needed_im_effective / CLOCKS_PER_SEC << '\n';
		cout << "\n\nNumber of matrix-vector multiplications:\n\n";
		cout << "Gradient descent: " << mult_count_gd << '\n';
		cout << "Iterative method: " << mult_count_im << '\n';
		cout << "Gradient descent (multiplication is implemented efficiently for a matrix of my type): " << mult_count_gd_effective << '\n';
		cout << "Iterative method (multiplication is implemented efficiently for a matrix of my type): " << mult_count_im_effective << '\n';
	}

	void diagram_data(int n)
	{
		My_type_matrix A(500);
		Gradient_descent<My_type_matrix> gd(A, true);
		Iterative_method<My_type_matrix> im(A, true);
	}

	void power(int n)
	{
		My_type_matrix A(n);
		Power_iter_max<My_type_matrix> pi_max(A);
		Power_iter_min pi_min(A);

		cout << "\n\n---------------------------Power iteration for n = " << n << "---------------------------------\n\n\n";

		cout << "Maximum eigenvalue (in absolute): " << pi_max.eigenvalue << '\n';
		cout << "Residual for the corresponding eigenvector: " << pi_max.norm_exp() << '\n';

		cout << "Minimum eigenvalue (in absolute): " << pi_min.eigenvalue << '\n';
		cout << "Residual for the corresponding eigenvector: " << pi_min.norm_exp() << '\n';
	}

	void newtons()
	{
		Newtons_method_for_my_function nm(1., 2., 1., 0.);
		nm.x.show();
	}

	void interpolation()
	{
		double a = -2, b = 2;

		ofstream fout("function.txt");
		fout << fixed << setprecision(10);
		double x = a, h_to_paint = (b - a) / 999;
		for (int i = 0; i < 1000; ++i)
		{
			func_values[i] = func(x);
			fout << x << '\t' << func_values[i] << '\n';
			x += h_to_paint;
		}
		fout.close();

		fout.open("errors.txt");
		ofstream fout_time("time.txt");

		double start_time;

		for (int i = 10; i <= 100; i += 10)
		{
			fout << i << '\t';
			fout_time << i << '\t';
			start_time = clock();
			Newton_equidistant ne(i, a, b, &func);
			fout_time << (clock() - start_time) / CLOCKS_PER_SEC * 1000 << '\t';
			fout << ne.max_d << '\t';
			start_time = clock();
			Newton_Chebyshev nc(i, a, b, &func);
			fout_time << (clock() - start_time) / CLOCKS_PER_SEC * 1000 << '\t';
			fout << nc.max_d << '\t';
			start_time = clock();
			Cubic_spline cs(i, a, b, &func);
			fout_time << (clock() - start_time) / CLOCKS_PER_SEC * 1000 << '\t';
			fout << cs.max_d;
			fout << '\n';
			fout_time << '\n';
		}
		fout.close();
		fout_time.close();
	}



	double a = -1, b = 1;
	double eps1 = 0.0001, eps2 = 0.000001, eps3 = 0.00000001;

	double exactValue = 0.436776185235761;


	void trapezoidal_rule()
	{
		cout << fixed << setprecision(12);

		cout << "Exact answer:\t" << exactValue << "\n\n";

		cout << "\nTrapezoidal rule:\n\n";
		numberF = 0;
		Trapezoidal_rule tr1(&f, &p, a, b, eps1);
		cout << "eps = " << tr1.eps << "\tanswer = " << tr1.ans << "\t\tn = " << tr1.n << "\t\tnumber of calculating f = " << numberF <<
			"\ndifference from exact value = " << abs(exactValue - tr1.ans) << "\n\n";
		numberF = 0;
		Trapezoidal_rule tr2(&f, &p, a, b, eps2);
		cout << "eps = " << tr2.eps << "\tanswer = " << tr2.ans << "\t\tn = " << tr2.n << "\tnumber of calculating f = " << numberF <<
			"\ndifference from exact value = " << abs(exactValue - tr2.ans) << "\n\n";
		numberF = 0;
		Trapezoidal_rule tr3(&f, &p, a, b, eps3);
		cout << "eps = " << tr3.eps << "\tanswer = " << tr3.ans << "\t\tn = " << tr3.n << "\tnumber of calculating f = " << numberF <<
			"\ndifference from exact value = " << abs(exactValue - tr3.ans) << "\n\n";
	}
};


int main()
{
	cout << fixed << setprecision(8);


	Program();



	cout << "\n\n\n";
	system("pause");
	return 0;
}