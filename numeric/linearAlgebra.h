#ifndef LINEARALGEBRA_H_
#define LINEARALGEBRA_H_

// my stuff
#include <other/Exception.h>


/**
 * Standard exception for this module.
 */
class LinearAlgebraException : public Exception {
public:
    LinearAlgebraException(const std::string & msg)
    	: Exception(msg)
    { }

    virtual ~LinearAlgebraException() throw()
    { }
};

template <typename Matrix, typename Vector>
Vector
powerIteration(const Matrix& A, const Vector& initVec, typename Vector::value_type tolerance = 1e-6, std::size_t maxIter = 30);


#endif /*LINEARALGEBRA_H_*/
