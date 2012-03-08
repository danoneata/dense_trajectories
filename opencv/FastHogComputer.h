#ifndef FASTHOGCOMPUTER_H_
#define FASTHOGCOMPUTER_H_

// STL
#include <vector>

// boost
#include <boost/numeric/ublas/vector.hpp>

// my stuff
#include <opencv/IplImageWrapper.h>
#include <opencv/IplImagePyramid.h>
#include <opencv/functions.h>
#include <geometry/Box.h>


class FastHogComputer
{
public:
	typedef std::vector<IplImagePyramid> IplImagePyramidVector;
	typedef double ValueType;
	typedef boost::numeric::ublas::vector<ValueType> VectorType;

protected:
	IplImagePyramidVector _integralHogBins;
	bool _fullOrientation;
	double _pyramidScaleFactor;

public:
	FastHogComputer(const IplImageWrapper& img, std::size_t nBins, bool fullOrientation = false, double pyramidScaleFactor = 0, bool imgIs2DVector = false);

	FastHogComputer(const IplImagePyramid& imgPyr, std::size_t nBins, bool fullOrientation = false);

	FastHogComputer(const FastHogComputer& fastHog);

	~FastHogComputer();

	FastHogComputer& operator=(const FastHogComputer& hogComputer);

	std::size_t numOfBins() const;

	std::size_t numOfLevels() const;

	double getPyramidScaleFactor() const;

	bool hasFullOrientation() const;

	VectorType getHog(const Box<int>& box, double scaleFactor) const
	{
		// get the pyramid level we are going to operate on
		std::size_t iLevel = _integralHogBins[0].getIndex(scaleFactor, -1);
		double scaleFactorLevel = _integralHogBins[0].getScaleFactorInv(iLevel);

		// scale the box according to the scale factors
		Box<int> scaledBox = scaleFactorLevel * box;
		return getHog(scaledBox, iLevel);
	}
	
	/**
	 * this is fast .. requires the coordinate to be adapted to the scale!
	 */
	VectorType getHog(const Box<int>& box_, std::size_t iLevel = 0) const;
	
	/**
	 * this is fast .. requires the coordinate to be adapted to the scale!
	 */
	ValueType getHogBin(const Box<int>& box, std::size_t iBin, std::size_t iLevel = 0) const
	{
		assert(iBin < _integralHogBins.size());
		return getIntegralRegion(_integralHogBins[iBin].getImage(iLevel), box);
	}

	/**
	 * this is fast .. requires the coordinate to be adapted to the scale!
	 */
	VectorType getIntegralHog(int x, int y, std::size_t iLevel = 0) const;
	
	std::size_t getSizeInBytes() const;
	
	const IplImagePyramidVector& getImagePyramidVector() const
	{
		return _integralHogBins;
	}

private:
	void init(const IplImageWrapper& orgImg, std::size_t nBins, bool fullOrientation, double pyramidScaleFactor, bool imgIs2DVector = false);

	void init(const IplImagePyramid& imgPyr, std::size_t nBins, bool fullOrientation, bool imgIs2DVector = false);
};

#include "FastHogComputer.hpp"

#endif /*FASTHOGCOMPUTER_H_*/
