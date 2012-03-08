#include "SumComputer.h"

// std
#include <iostream>
#include <exception>

// boost
#include <boost/lexical_cast.hpp>

// opencv
#include <opencv/cv.h>

// my stuff
#include <opencv/functions.h>


SumComputer::SumComputer(const IplImageWrapper& orgImg, double pyramidScaleFactor)
{
	if (!(orgImg->nChannels == 1 || orgImg->nChannels == 3))
		throw std::runtime_error("SumComputer(): the image has to be either a gray or a color image!");
	init(orgImg, pyramidScaleFactor);
}

SumComputer::SumComputer(const IplImagePyramid& imgPyr)
{
	if (imgPyr.getImage((std::size_t)0)->nChannels != 1)
		throw std::runtime_error("SumComputer(): imgPyr is not a gray image!");
	if (imgPyr.getImage((std::size_t)0)->depth != IPL_DEPTH_64F)
		throw std::runtime_error("SumComputer(): format of imgPyr is not double!");
	init(imgPyr);
}

SumComputer::~SumComputer()
{

}

void SumComputer::init(const IplImageWrapper& orgImg, double pyramidScaleFactor)
{
	// convert image to gray float image
	IplImageWrapper img = convert2GrayFloatImg(orgImg);
	
	// build up a scale pyramid
	IplImagePyramid imgPyr(img, pyramidScaleFactor);
	img = IplImageWrapper();

	// continue in other init method
	init(imgPyr);
}

void SumComputer::init(const IplImagePyramid& imgPyr)
{
	// create a pyramid for of integral images
	_integralImgPyr = IplImagePyramid(cvGetSize(imgPyr.getImage((std::size_t)0)), IPL_DEPTH_64F, 1, imgPyr.getScaleFactor());
	assert(_integralImgPyr.numOfLevels() == imgPyr.numOfLevels());

	// compute for each image in the pyramid its integral representation
	for (std::size_t i = 0; i < imgPyr.numOfLevels(); ++i) {
		IplImageWrapper img = imgPyr.getImage(i);
		
		// if the image has a mask, we will set all pixel values
		// outside the mask to zero
		if (img.hasMask()) {
			// make a deep copy of the image at the current level and set 
			// all values outside the mask to zero
			img = img.clone();
			setZeroOutsideMask(img, img.getMask());
		}
		computeIntegralImage(img, _integralImgPyr.getImage(i));
	}
}

std::size_t SumComputer::getSizeInBytes() const
{
	std::size_t size(0);
	for (std::size_t j = 0; j < _integralImgPyr.numOfLevels(); ++j)
		size += _integralImgPyr.getImage(j)->height * _integralImgPyr.getImage(j)->widthStep;
	return size;
}

