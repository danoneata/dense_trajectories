#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

#include <other/Exception.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

/**
 * Standard exception for this module.
 */
class OptimizationException : public Exception {
public:
    OptimizationException(const std::string & msg)
    	: Exception(msg)
    { }

    virtual ~OptimizationException() throw()
    { }
};

/**
 * Given a starting point p[1..n] that is a vector of length n, the Broyden-Fletcher-Goldfarb-
 * Shanno variant of Davidon-Fletcher-Powell minimization is performed on a function func, using
 * its gradient as calculated by a routine dfunc. The convergence requirement on zeroing the
 * gradient is input as gtol. Returned quantities are p[1..n] (the location of the minimum),
 * iter (the number of iterations that were performed), and fret (the minimum value of the
 * function). The routine lnsrch is called to perform approximate line minimizations.
 *
 * func: takes as parameter one vector where it evaluates the vector's value
 * dfunc: takes as parameter one vector and returns the gradient at the given vector
 *
 * OptimizationFunction is an functor evaluating a new position with the following form:
 *
 * struct OptimizationFunctor {
 *     double operator()(const Vector& newPos);
 * }
 *
 * GradientFunction is an functor computing a gradient vector at a given point,
 * it has the following form:
 *
 * struct GradientFunctor {
 *     Vector operator()(const Vector& pos);
 * }
 *
 * (adapted from numerical recipies in C)
 */
template<typename Vector, typename OptimizationFunction, typename GradientFunction>
void minimizationBFGS(Vector& p, typename Vector::value_type gtol, int *iter, typename Vector::value_type *fret,
		OptimizationFunction& func, GradientFunction& dfunc)
		throw (OptimizationException);

/**
 * Given an n-dimensional point xold[1..n], the value of the function and gradient there, fold
 * and g[1..n], and a direction p[1..n], finds a new point x[1..n] along the direction p from
 * xold where the function func has decreased "sufficiently." The new function value is returned
 * in f. stpmax is an input quantity that limits the length of the steps so that you do not try to
 * evaluate the function in regions where it is undefined or subject to overflow. p is usually the
 * Newton direction. The output quantity check is false (0) on a normal exit. It is true (1) when
 * x is too close to xold. In a minimization algorithm, this usually signals convergence and can
 * be ignored. However, in a zero-finding algorithm the calling program should check whether the
 * convergence is spurious. Some "difficult" problems may require double precision in this routine.
 *
 * func: takes as parameter one vector where it evaluates the vector's value
 *
 * (adapted from numerical recipies in C)
 */
template<typename Vector, typename OptimizationFunction>
void lineSearch(const Vector& xold, typename Vector::value_type fold, const Vector& g, Vector& p, Vector& x,
		typename Vector::value_type* f, typename Vector::value_type stpmax,
		bool* ok, OptimizationFunction& func)
		throw (OptimizationException);


#include "optimization.hpp"

#endif // OPTIMIZATION_H
