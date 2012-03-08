#include "clustering.h"

// std libs
#include <algorithm>
#include <list>

// boost libs
#include <boost/typeof/typeof.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

// my libs
#include <boost/ublasFunctions.hpp>
//#include <other/functions.h>

// namespaces
namespace ublas = boost::numeric::ublas;
using std::vector;
using std::list;
using std::cout;
using std::cerr;
using std::endl;

std::vector<boost::numeric::ublas::vector<double> >
meanShiftClustering(const std::vector<boost::numeric::ublas::vector<double> >& dataPoints, const std::vector<boost::numeric::ublas::vector<double> >& variances,
		const boost::numeric::ublas::vector<double>& weights, const boost::numeric::ublas::vector<double>& stopValues,
		boost::numeric::ublas::vector<double>* _modeWeights, boost::numeric::ublas::vector<std::size_t>* _mapping)
{
	assert(dataPoints.size() == variances.size());
	assert(dataPoints.size() == weights.size());

	// if we have not data points it's easy :)
	if (dataPoints.size() <= 0) {
		if (_modeWeights)
			_modeWeights->resize(0u);
		if (_mapping)
			_mapping->resize(0u);
		return vector<boost::numeric::ublas::vector<double> >(0);
	}

	// get the dimensionality and the size of the data set
	std::size_t nDims = dataPoints.front().size();
	std::size_t nPoints = dataPoints.size();
	assert(stopValues.size() == nDims);

	// create a list of uncertainty matrices, for each data point one
	vector<ublas::mapped_matrix<double> > covMatrices(nPoints);
	vector<ublas::mapped_matrix<double> > covMatricesInv(nPoints);
	ublas::vector<double> constSqrtDets(nPoints);

	// compute the uncertainty matrices, their inverse for all detections,
	// and a particular constant
	for (std::size_t i = 0; i < nPoints; ++i) {
		assert(dataPoints[i].size() == nDims);
		assert(variances[i].size() == nDims);

		// initiate the uncertainty matrix and its invers (can be computed directly
		// since we have a diagonal matrix)
		covMatrices[i] = ublas::mapped_matrix<double>(nDims, nDims, nDims);
		covMatricesInv[i] = ublas::mapped_matrix<double>(nDims, nDims, nDims);
		for (std::size_t iDim = 0; iDim < nDims; ++iDim) {
			covMatrices[i](iDim, iDim) = variances[i][iDim];
			covMatricesInv[i](iDim, iDim) = 1 / variances[i][iDim];
		}

		// compute the constant: 1 / sqrt(det(M))
		constSqrtDets[i] = 1 / sqrt(matrixDeterminant(covMatrices[i]));
	}

	// compute the mode (i.e., its maximum) for each detection
	ublas::vector<double> tmpWeights(nPoints);
	ublas::vector<double> mode(nDims);
	ublas::vector<double> newMode(nDims);
	ublas::vector<double> tmpDiff(nDims);
	ublas::mapped_matrix<double> meanCovMatrix(nDims, nDims, nDims);
	ublas::mapped_matrix<double> meanCovMatrixInv(nDims, nDims, nDims);
	vector<ublas::vector<double> > tmpModes(nPoints);
	ublas::vector<double> tmpModeWeights(nPoints);
	vector<ublas::mapped_matrix<double> > tmpCovMatricesInv(nPoints);
	ublas::mapped_matrix<double> tmpMatrix(nDims, nDims, nDims);
	ublas::vector<double> tmpVec(nDims);
//cout << "### new clustering ###" << endl;
	for (std::size_t i = 0; i < nPoints; ++i) {
//cout << "=== new Point ===" << endl;
//cout << "dataPoint: " << dataPoints[i] << endl;
		// compute the current mode iteratively
		mode.assign(dataPoints[i]);
//cout << "mode: " << mode << endl;
//list<ublas::vector<double> > tmpVecList;
//tmpVecList.push_back(mode);
//list<ublas::matrix<double> > tmpMatList;
		double weight = weights[i];
		bool isDone(false);
		while (!isDone) {
			// compute the weight of each detection given the current mode
			for (std::size_t j = 0; j < nPoints; ++j) {
				ublas::noalias(tmpDiff) = mode - dataPoints[j];
				ublas::axpy_prod(covMatricesInv[j], tmpDiff, tmpVec, true);
				tmpWeights[j] = constSqrtDets[j] * weights[j]
						* exp(-0.5 * ublas::inner_prod(tmpDiff, tmpVec));
			}

			// normalize the weights such that their sum equals 1
			ublas::vector<double> oldWeights = tmpWeights;
			double sumTmpWeights = ublas::norm_1(tmpWeights);
			if (sumTmpWeights == 0.0) {
				isDone = true;
				continue;
			}
			tmpWeights *= 1 / sumTmpWeights;
//if (isnan(ublas::norm_1(tmpWeights))) {
//	cout << "ERROR NAN A2" << endl;
//	cout << "    " << weights << endl;
//}

			// compute the weighted harmonic mean convariance matrix
			// (since we have diagonal matrices, the inverse can be computed directly)
			meanCovMatrixInv.clear();
			meanCovMatrix.clear();
			for (std::size_t j = 0; j < nPoints; ++j)
				ublas::noalias(meanCovMatrixInv) += tmpWeights[j] * covMatricesInv[j];
			for (ublas::mapped_matrix<double>::const_iterator1 i1 = meanCovMatrixInv.begin1(); i1 != meanCovMatrixInv.end1(); ++i1)
				for (ublas::mapped_matrix<double>::const_iterator2 i2 = i1.begin(); i2 != i1.end(); ++i2)
					meanCovMatrix(i2.index1(), i2.index2()) = 1 / *i2;
//tmpMatList.push_back(meanCovMatrixInv);
//tmpMatList.push_back(meanCovMatrix);

			// update the position of the mode
			newMode.assign(ublas::zero_vector<double>(nDims));
			for (std::size_t j = 0; j < nPoints; ++j) {
				ublas::axpy_prod(covMatricesInv[j], dataPoints[j], tmpVec, true);
				ublas::noalias(newMode) += tmpWeights[j] * tmpVec;
			}
			newMode = prod(meanCovMatrix, newMode);
//tmpVecList.push_back(tmpWeights);
//tmpVecList.push_back(newMode);

			// compute the distance between the last mode position and the new one
			// and check whether we are done are not
			ublas::noalias(tmpDiff) = newMode - mode;
			isDone = true;
			for (std::size_t iDim = 0; iDim < nDims && isDone; ++iDim)
				if (fabs(tmpDiff[iDim]) > fabs(stopValues[iDim]))
					isDone = false;

			// update the mode position
			mode.assign(newMode);
//tmpVecList.push_back(mode);
		};
		
//for(BOOST_AUTO(iVec, tmpVecList.begin()); iVec != tmpVecList.end(); ++iVec)
//	cout << *iVec << endl;
//cout << endl;
//
//for(BOOST_AUTO(iMat, tmpMatList.begin()); iMat != tmpMatList.end(); ++iMat)
//	cout << *iMat << endl << endl;
//cout << endl;

		// done .. we found the mode for the current point, save it
		tmpModes[i] = mode;
		tmpModeWeights[i] = weight;
		tmpCovMatricesInv[i] = meanCovMatrixInv;
	}

	//
	// fuze tmpModes together
	//

	ublas::vector<std::size_t> mapping(nPoints);
	std::fill(mapping.begin(), mapping.end(), nPoints);
	ublas::vector<double> distances = ublas::zero_vector<double>(nPoints);

	// go through the list and label tmpModes that need to be fuzed
	std::size_t nModes(0u);
	for (std::size_t i = 0; i < nPoints; ++i) {

		// if the current mode does not have a label, we found a new mode
		// otherwise it already has a label and has already been assigned
		if (mapping[i] >= nPoints) {
			mapping[i] = nModes;
			++nModes;
		}

		// go through the list and merge current averaged mode with other tmpModes
		for (std::size_t j = 0; j < nPoints; ++j) {
			if (i == j)
				continue;

			// compute the Mahalanobis distance between the two tmpModes
			ublas::noalias(tmpDiff) = tmpModes[i] - tmpModes[j];
			ublas::axpy_prod(tmpCovMatricesInv[i], tmpDiff, tmpVec, true);
			double tmpDistance = inner_prod(tmpDiff, tmpVec);
			if (tmpDistance < nDims && (mapping[j] >= nPoints || tmpDistance < distances[j])) {
				mapping[j] = mapping[i];
				distances[j] = tmpDistance;
			}
		}
	}

//cout << "###### TMP MODES #######" << endl;
//for (std::size_t i = 0; i < nPoints; ++i)
//	cout << mapping[i] << " >> " << distances[i] << " -- " << tmpModes[i] << endl;

	// compute the real modes
	ublas::vector<double> modeWeights(nModes);
	ublas::vector<double> sumWeights(nModes);
	ublas::vector<std::size_t> nPointsPerMode(nModes);
	sumWeights.assign(ublas::zero_vector<double>(nModes));
	modeWeights.assign(ublas::zero_vector<double>(nModes));
//	nPointsPerMode.assign(ublas::zero_vector<std::size_t>(nModes));
	std::fill(nPointsPerMode.begin(), nPointsPerMode.end(), 0);
	std::vector<boost::numeric::ublas::vector<double> > modes(nModes);
	for (std::size_t i = 0; i < nModes; ++i)
		modes[i] = ublas::zero_vector<double>(nDims);

	// compute for each mode the sum of weights of all points that belong to it
	for (std::size_t i = 0; i < nPoints; ++i) {
		std::size_t iMode = mapping[i];
		assert(nModes > iMode);

		modeWeights[iMode] = std::max(modeWeights[iMode], tmpModeWeights[i]);
//		modeWeights[iMode] += tmpModeWeights[i] * tmpModeWeights[i];
		sumWeights[iMode] += tmpModeWeights[i];
		++nPointsPerMode[iMode];
	}

	// collect the points
	for (std::size_t i = 0; i < nPoints; ++i) {
		std::size_t iMode = mapping[i];
		assert(nModes > iMode);

		// if all points have the weight 0.0, compute the average mode position
		// otherwise the weighted mean mode position
		if (0.0 == sumWeights[iMode])
			ublas::noalias(modes[iMode]) += tmpModes[i];
		else
			ublas::noalias(modes[iMode]) += tmpModeWeights[i] * tmpModes[i];
	}

	// compute the mean mode
	for (std::size_t i = 0; i < nModes; ++i) {
		// if all points have the weight 0.0, compute the average mode position
		// otherwise the weighted mean mode position
		if (0.0 == sumWeights[i])
			modes[i] *= 1.0 / nPointsPerMode[i];
		else
			modes[i] *= 1.0 / sumWeights[i];
//		modeWeights[i] *= 1 / sumWeights[i];
	}

	// return the results
	if (_modeWeights)
		*_modeWeights = modeWeights;
	if (_mapping)
		*_mapping = mapping;
	return modes;
}




