// standard stuff
#include <cmath>
#include <iostream>
#include <cassert>
//#include <boost/numeric/ublas/io.hpp>

template<typename Vector, typename OptimizationFunction>
void lineSearch(const Vector& xold, typename Vector::value_type fold, const Vector& g, Vector& p, Vector& x,
		typename Vector::value_type* f, typename Vector::value_type stpmax,
		bool* ok, OptimizationFunction& func)
		throw (OptimizationException)
{
	using namespace boost::numeric::ublas;
	typedef typename Vector::value_type value_type;

	assert(norm_1(p) > 0);
	assert(norm_1(g) > 0);
	assert(stpmax > 0);

	// Ensures sufficient decrease in function value.
	value_type ALF = 1.0e-4;
	// Convergence criterion on delta_x.
	value_type TOLX = 1.0e-7;
	int n = xold.size();
	int i;
	value_type a,alam,alam2(1.0),alamin,b,disc,f2(*f),rhs1,rhs2,slope,sum,temp,test,tmplam;
	*ok = true;

	// Scale if attempted step is too big.
	sum = norm_2(p);
	if (sum > stpmax)
		p *= stpmax/sum;

	// check whether gradient and search direction are orthogonal
	slope = inner_prod(g, p);
	if (slope > 0.0)
		throw OptimizationException("Roundoff problem in lineSearch.");

	// Compute lamda_min
	test = 0.0;
	for (i = 0; i < n; ++i) {
		temp = fabs(p[i]) / std::max<value_type>(fabs(xold[i]), 1.0);
		if (temp > test)
			test = temp;
	}
	alamin = TOLX / test;

	// Always try full Newton step first. Start of iteration loop.
	alam = 1.0;
	for (;;) {
		x = xold + alam * p;
		*f = func(x);

		if (alam < alamin) {
			// Convergence on delta-x. For zero finding, the calling program should
			// verify the convergence.
			x = xold;
			*ok = false;
			return;
		}
		else if (*f <= fold + ALF * alam * slope) {
			// Sufficient function decrease.
			return;
		}
		else if (std::isinf(*f) || std::isnan(*f)) {
			// we went too far .. choose a smaller alam
			alam = alam * 0.1;
			continue;
		}
		else {
			// Backtrack.
			if (alam == 1.0)
				// First time.
				tmplam = -slope / (2.0 * (*f - fold - slope));
			else {
				// Subsequent backtracks.
				rhs1 = *f - fold - alam * slope;
				rhs2 = f2 - fold - alam2 * slope;
				a = (rhs1 / (alam * alam) - rhs2 / (alam2 * alam2)) / (alam - alam2);
				b = (-alam2 * rhs1 / (alam * alam) + alam * rhs2 / (alam2 * alam2)) / (alam - alam2);
				if (a == 0.0)
					tmplam = -slope / (2.0 * b);
				else {
					disc = b * b -3.0 * a * slope;
					if (disc < 0.0)
						tmplam = 0.5 * alam;
					else if (b <= 0.0)
						tmplam = (-b + sqrt(disc)) / (3.0 * a);
					else
						tmplam = -slope / (b + sqrt(disc));
				}
				if (tmplam > 0.5 * alam)
					tmplam = 0.5 * alam;
			}
		}
		alam2 = alam;
		f2 = *f;
		alam = std::max<value_type>(tmplam, 0.1 * alam);
		// Try again.
	}
}

template<typename Vector, typename OptimizationFunction, typename GradientFunction>
void minimizationBFGS(Vector& p, typename Vector::value_type gtol, int *iter, typename Vector::value_type *fret,
		OptimizationFunction& func, GradientFunction& dfunc)
		throw (OptimizationException)
{
	using namespace boost::numeric;
	typedef typename Vector::value_type value_type;

	// Maximum allowed number of iterations.
	const value_type ITMAX = 200;
	// Machine precision.
	const value_type EPS = 3.0e-8;
	// Convergence criterion on x values.
	const value_type TOLX = 4 * EPS;
	// Scaled maximum step length allowed in line searches.
	const value_type STPMX = 100.0;

	bool ok;
	int i, its, j;
	value_type den, fac, fad, fae, fp, stpmax, sumdg, sumxi, temp, test;
	int n = p.size();
	ublas::vector<value_type> dg(n), g(n), hdg(n);
	ublas::vector<value_type> pnew(n), vector(n), xi(n);
	ublas::matrix<value_type> hessin(n, n);

	hessin.assign(ublas::identity_matrix<value_type>(n));
//std::cout << "#OPT p: " << p << std::endl;
	g = dfunc(p);
//std::cout << "#OPT g: " << g << std::endl;
	xi = -g;
	fp = func(p);
//std::cout << "#OPT fp: " << fp << std::endl;
	stpmax = STPMX * std::max<value_type>(ublas::norm_2(p), n);
	*iter = 0;
	*fret = fp;

	// in case the computed gradient vector is a zero-vector, we exit the function
	if (ublas::norm_1(g) == 0 || stpmax <= 0)
		return;

	// Main loop over the iterations.
	for (its = 0; its < ITMAX; ++its) {
		*iter = its + 1;
		lineSearch(p, fp, g, xi, pnew, fret, stpmax, &ok, func);

		// The new function evaluation occurs in lineSearch; save the function value in fp for the
		// next line search. It is usually safe to ignore the value of "ok".
		fp = *fret;

		// Update the line direction, and the current point.
		xi = pnew - p;
		p = pnew;

		// Test for convergence on x.
		test = 0.0;
		for (i = 0; i < n; ++i) {
			temp = fabs(xi[i]) / std::max<value_type>(fabs(p[i]), 1.0);
			if (temp > test)
				test = temp;
		}
		if (test < TOLX) {
			return;
		}

		// Save the old gradient, and get the new gradient.
		dg = g;
		g = dfunc(p);

		// Test for convergence on zero gradient.
		test = 0.0;
		den = std::max<value_type>(*fret, 1.0);
		for (i = 0; i < n; ++i) {
			temp = fabs(g[i]) * std::max<value_type>(fabs(p[i]), 1.0) / den;
			if (temp > test)
				test=temp;
		}
		if (test < gtol) {
			return;
		}

		// Compute difference of gradients, and difference times current matrix.
		dg = g - dg;
		hdg = ublas::prod(hessin, dg);

		// Calculate dot products for the denominators.
		fac = ublas::inner_prod(dg, xi);
		fae = ublas::inner_prod(dg, hdg);
		sumdg = ublas::inner_prod(dg, dg);
		sumxi = ublas::inner_prod(xi, xi);

		// Skip update if fac not suffciently posifactive.
		if (fac > sqrt(EPS * sumdg * sumxi)) {
			fac = 1.0 / fac;
			fad = 1.0 / fae;

			// The vector that makes BFGS different from DFP:
			dg = fac * xi - fad * hdg;

			// The BFGS updating formula:
			for (i = 0; i < n; i++) {
				for (j = i; j < n; j++) {
					hessin(i, j) += fac * xi[i] * xi[j] - fad * hdg[i] * hdg[j] + fae * dg[i] * dg[j];
					hessin(j, i) = hessin(i, j);
				}
			}
		}
		// Now calculate the next direction to go, and go back for another iteration.
		xi = -ublas::prod(hessin, g);
	}
	throw OptimizationException("too many iterations in dfpmin");
}
