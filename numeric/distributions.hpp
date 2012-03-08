#ifndef DISTRIBUTIONS_HPP_
#define DISTRIBUTIONS_HPP_

#include "distributions.h"

// std libs
#include <algorithm>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

// namespaces
namespace l = boost::lambda;


template <class InputIterator>
Gaussian Gaussian::create(InputIterator first, InputIterator last, bool isNormalized)
{
	std::size_t n = std::distance(first, last);
	assert(n > 1);

	double mean(0);
	std::for_each(first, last, mean += l::_1);
	mean /= n;

	double variance(0);
	std::for_each(first, last, variance += l::bind(&pow, l::_1 - mean, 2));
	variance /= (n - 1);

	double sigma = sqrt(variance);

	return Gaussian(mean, sigma, isNormalized);
}

#endif
