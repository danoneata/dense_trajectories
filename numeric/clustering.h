#ifndef CLUSTERING_H_
#define CLUSTERING_H_

#include <vector>
#include <boost/numeric/ublas/vector.hpp>

/**
 * Performs on a given set of data points of dimensionality D a Mean-Shift clustering
 * with variable band-width.
 *
 * @param dataPoints  N data points of dimensionality D
 * @param variances   gives variances as uncertainty for each dimension and for each
 *                    data point (note that the uncertainties must be non-zero!)
 * @param weights     weights for each data point
 * @param stopValues  computation of a mode is stopped as the difference between the current
 *                    and the last position is in all dimensions smaller than this vector
 * @param weights     if not NULL, returns a vector with the weights for each mode
 * @param mapping     if not NULL, returns a vector that indicates the mapping/assignment
 *                    of data point to cluster
 * @return            a vector of mode centers
 */
std::vector<boost::numeric::ublas::vector<double> >
meanShiftClustering(const std::vector<boost::numeric::ublas::vector<double> >& dataPoints, const std::vector<boost::numeric::ublas::vector<double> >& variances,
		const boost::numeric::ublas::vector<double>& weights, const boost::numeric::ublas::vector<double>& stopValues,
		boost::numeric::ublas::vector<double>* modeWeights = 0, boost::numeric::ublas::vector<std::size_t>* modeMapping = 0);

#endif /*CLUSTERING_H_*/
