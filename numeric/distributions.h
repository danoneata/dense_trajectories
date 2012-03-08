#ifndef DISTRIBUTIONS_H_
#define DISTRIBUTIONS_H_

// std libs
#include <cmath>

class Gaussian {
protected:
	double _mu;
	double _sigma;
	bool _isNormalized;
	mutable double _2sigmaSquareInv;
	mutable double _const;
	mutable bool _areValuesCached;

public:
	Gaussian()
		: _mu(0), _sigma(1), _isNormalized(false),
		_2sigmaSquareInv(0), _const(1), _areValuesCached(true)
	{ }

	// Gaussian(const Gaussian& gaussian)

	Gaussian(double mu, double sigma, bool isNormalized = true)
		: _mu(mu), _sigma(sigma), _isNormalized(isNormalized),
		_2sigmaSquareInv(0), _const(1), _areValuesCached(false)
	{ }

	// Gaussian& operator=(const Gaussian& gaussian)

	double operator()(double x) const
	{
		// check whether we already computed the cached values
		if (!_areValuesCached)
			cacheValues();

		// compute the value for the Gaussian distribution
		double xMinMu = x - _mu;
		double val = std::exp(-xMinMu * xMinMu * _2sigmaSquareInv);
		if (_isNormalized)
			return _const * val;
		return val;
	}

	double getMu() const
	{
		return _mu;
	}

	void setMu(double mu)
	{
		_mu = mu;
		_areValuesCached = false;
	}

	double getSigma() const
	{
		return _sigma;
	}

	void setSigma(double sigma)
	{
		_sigma = sigma;
		_areValuesCached = false;
	}

	bool isNormalized() const
	{
		return _isNormalized;
	}

	void setNormalized(bool isNormalized)
	{
		_isNormalized = isNormalized;
	}

	template <class InputIterator>
	static
	Gaussian create(InputIterator first, InputIterator last, bool isNormalized = true);

protected:
	void cacheValues() const
	{
		_2sigmaSquareInv = 1 / (2 * _sigma * _sigma);
		_const = 1 / (_sigma * sqrt(2 * M_PI));
		_areValuesCached = true;
	}

};

#endif /*DISTRIBUTIONS_H_*/
