#ifndef SUMCOMPUTER_H_
#define SUMCOMPUTER_H_

// my stuff
#include <opencv/IplImageWrapper.h>
#include <opencv/IplImagePyramid.h>
#include <opencv/functions.h>
#include <geometry/Box.h>


class SumComputer
{
protected:
	IplImagePyramid _integralImgPyr;

public:
	SumComputer(const IplImageWrapper& img, double pyramidScaleFactor = 0);

	SumComputer(const IplImagePyramid& imgPyr);

	~SumComputer();

	double getSum(const Box<int>& box, double scaleFactor) const
	{
		// get the scale Level and scale the box correctly
		std::size_t iLevel = _integralImgPyr.getIndex(scaleFactor, -1);
		Box<int> scaledBox(box);
		scaledBox.scale(_integralImgPyr.getScaleFactorInv(iLevel), _integralImgPyr.getScaleFactorInv(iLevel));
		
		// pass to the other getSum() method
		return getSum(scaledBox, iLevel);
	}
	
	/**
	 * this is fast .. requires the coordinate to be adapted to the scale!
	 */
	double getSum(const Box<int>& box, std::size_t iLevel = 0) const
	{
		return getIntegralRegion(_integralImgPyr.getImage(iLevel), box);
	}
	
	std::size_t getSizeInBytes() const;
	
	const IplImagePyramid& getImagePyramid() const
	{
		return _integralImgPyr;
	}

private:
	void init(const IplImageWrapper& img, double pyramidScaleFactor);

	void init(const IplImagePyramid& imgPyr);
};

#endif /*SUMCOMPUTER_H_*/
