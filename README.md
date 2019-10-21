# NumericalAnalysis
Asembly of numerical methods implemented in C++.

Solving a linear system: [Gaussian elimination](https://github.com/yanava99/NumericalAnalysis/blob/master/Gaussian_elimination.h), [iterative method](https://github.com/yanava99/NumericalAnalysis/blob/master/Iterative_method.h), [Newton's method](https://github.com/yanava99/NumericalAnalysis/blob/master/Newtons_method.h).

Interpolation: [Cubic spline](https://github.com/yanava99/NumericalAnalysis/blob/master/Cubic_spline.h), [Newton polynomial](https://github.com/yanava99/NumericalAnalysis/blob/master/Newton_polynomial.h) parent class with [Chebyshev nodes](https://github.com/yanava99/NumericalAnalysis/blob/master/Newton_polynomial_Chebyshev_nodes.h) and [equidistant nodes](https://github.com/yanava99/NumericalAnalysis/blob/master/Newton_polynomial_equidistant_nodes.h) inheritors.

Finding min of a function: [Gradient descent](https://github.com/yanava99/NumericalAnalysis/blob/master/Gradient_descent.h).

Finding min/max eigenvalue and corresponding eigenvector: [Power iteration](https://github.com/yanava99/NumericalAnalysis/blob/master/Power_iteration.h).

Approximating the definite integral: [Trapezoidal rule](https://github.com/yanava99/NumericalAnalysis/blob/master/Trapezoidal_rule.h).

[Vector](https://github.com/yanava99/NumericalAnalysis/blob/master/Vector.h) and [Matrix](https://github.com/yanava99/NumericalAnalysis/blob/master/Matrix.h) are supporting classes; [main](https://github.com/yanava99/NumericalAnalysis/blob/master/main.cpp) is used for testing.
