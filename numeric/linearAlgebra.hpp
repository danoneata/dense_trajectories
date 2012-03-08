#ifndef LINEARALGEBRA_HPP_
#define LINEARALGEBRA_HPP_

#include "linearAlgebra.h"

// boost stuff
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <iostream>
#include <boost/numeric/ublas/io.hpp>


template <typename Matrix, typename Vector>
Vector
powerIteration(const Matrix& A, const Vector& initVec, typename Vector::value_type tolerance, std::size_t maxIter)
{
	namespace ublas = boost::numeric::ublas;
//	namespace l = boost::lambda;
	using std::cout;
	using std::endl;
	typedef typename Vector::value_type ValueType;
//	typedef typename std::iterator_traits<FeatureVecIterator>::value_type WeightType;
//	typedef typename std::iterator_traits<FeatureVecIterator>::value_type VectorType;
//	typedef typename VectoryType::value_type ValueType;

	// check the norm of the starting vector
	ValueType initNorm = ublas::norm_2(initVec);
	if (initNorm < 1e-10)
		throw LinearAlgebraException("norm of the init vector is too small");

	Vector eigenVec(initVec.size());
	eigenVec = initVec;
	eigenVec *= 1 / initNorm;
	Vector newEigenVec(initVec.size());

//std::cout << "init: " << eigenVec << std::endl << std::endl;

	// iterate until the maximum of iterations is achieved or the solution converged
	for (std::size_t iIter = 0; iIter < maxIter; ++iIter) {
		newEigenVec = ublas::prod(A, eigenVec);
		newEigenVec *= 1 / ublas::norm_2(newEigenVec);
//std::cout << iIter << ") " << newEigenVec << std::endl;
//std::cout << "   " << eigenVec << std::endl;
//std::cout << "   " << ublas::norm_2(eigenVec - newEigenVec) << std::endl << std::endl;

		// check whether the solution converged
		if (ublas::norm_2(eigenVec - newEigenVec) < tolerance)
			break;

		eigenVec = newEigenVec;
	}

	return newEigenVec;
}

#endif
