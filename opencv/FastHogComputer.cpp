#include "FastHogComputer.h"

// std and STL
#include <cmath>
#include <fstream>
#include <iostream>

// opencv
#include <opencv/cv.h>

// my stuff
#include <opencv/functions.h>
#include <geometry/Box.hpp>

// tmp
#include <opencv/highgui.h>
#include <boost/lexical_cast.hpp>
using boost::lexical_cast;
using std::string;

using std::cerr;
using std::cout;
using std::endl;


void FastHogComputer::init(const IplImageWrapper& orgImg, std::size_t nBins, bool fullOrientation, double pyramidScaleFactor, bool imgIs2DVector)
{
	// convert image to gray float image
	IplImageWrapper img;
	if (imgIs2DVector)
        img = orgImg.clone();
    else
		img = convert2GrayFloatImg(orgImg);

	// make sure that 2D vectorial data is float data with two channels
	if (imgIs2DVector && (orgImg->depth != IPL_DEPTH_32F || orgImg->nChannels != 2))
		throw std::runtime_error("FastHogComputer::init() - 2D vectorial data is not in correct format");
	
	// build up a scale pyramid
	IplImagePyramid imgPyr(img, pyramidScaleFactor);
	img = IplImageWrapper();

	// continue in other init method
	init(imgPyr, nBins, fullOrientation, imgIs2DVector);
}


void FastHogComputer::init(const IplImagePyramid& imgPyr, std::size_t nBins, bool fullOrientation, bool imgIs2DVector)
{
	int nChannels = 1;
	CvSize size = cvGetSize(imgPyr.getImage((std::size_t)0));

	IplImagePyramid gradientOrientation(size, IPL_DEPTH_32F, 1, imgPyr.getScaleFactor());
	IplImagePyramid gradientMagnitude(size, IPL_DEPTH_32F, 1, imgPyr.getScaleFactor());
	IplImagePyramidVector hogBins(nBins);
	for (std::size_t i = 0; i < nBins; ++i) {
		_integralHogBins[i] = IplImagePyramid(size, IPL_DEPTH_64F, 1, imgPyr.getScaleFactor());
		hogBins[i] = IplImagePyramid(size, IPL_DEPTH_64F, 1, imgPyr.getScaleFactor());
	}

	// some other variables that we need
	float maxDeg = fullOrientation ? 360 : 180;
	float invBinWidthDeg = nBins / maxDeg;
	std::size_t nLevels = gradientOrientation.numOfLevels();
	assert (gradientOrientation.numOfLevels() == nLevels);
	assert (gradientMagnitude.numOfLevels() == nLevels);

	// iterate over the different pyramid levels
	for (std::size_t iLevel = 0; iLevel < nLevels; ++iLevel) {
		// get references to the current level in the pyramid
		IplImageWrapper gradientOrientationTmp = gradientOrientation.getImage(iLevel);
		IplImageWrapper gradientMagnitudeTmp = gradientMagnitude.getImage(iLevel);
		IplImageWrapper xDerivTmp(cvGetSize(gradientOrientationTmp), IPL_DEPTH_32F, nChannels);
		IplImageWrapper yDerivTmp(cvGetSize(gradientOrientationTmp), IPL_DEPTH_32F, nChannels);
		if (!imgIs2DVector) {
			cvSetZero(xDerivTmp);
			cvSetZero(yDerivTmp);
			cvSobel(imgPyr.getImage(iLevel), xDerivTmp, 1, 0, 1);
			cvSobel(imgPyr.getImage(iLevel), yDerivTmp, 0, 1, 1);
		}
		else {
			cvSplit(imgPyr.getImage(iLevel), xDerivTmp, yDerivTmp, NULL, NULL);
		}

		// compute orientation + magnitude
		char* rowXDeriv = xDerivTmp->imageData;
		char* rowYDeriv = yDerivTmp->imageData;
		char* rowOrientation = gradientOrientationTmp->imageData;
		char* rowMagnitude = gradientMagnitudeTmp->imageData;

		// iterate over rows
		int levelWidth = xDerivTmp->width;
		int levelHeight = xDerivTmp->height;
		for (int iy = 0; iy < levelHeight;
				++iy, rowXDeriv += xDerivTmp->widthStep, rowYDeriv += yDerivTmp->widthStep,
				rowOrientation += gradientOrientationTmp->widthStep, rowMagnitude += gradientMagnitudeTmp->widthStep) {
			float* xXDeriv = (float*) rowXDeriv;
			float* xYDeriv = (float*) rowYDeriv;
			float* xOrientation = (float*) rowOrientation;
			float* xMagnitude = (float*) rowMagnitude;

			// iterate over columns
			for (int ix = 0; ix < levelWidth;
					++ix, ++xOrientation, ++xMagnitude) {
				// iterate over channels
				// TODO: in order to deal with color information it should not be done like this
				for (int ic = 0; ic < nChannels; ++ic, ++xXDeriv, ++xYDeriv) {
					float magnitude = sqrtf(powf(*xYDeriv, 2) + powf(*xXDeriv, 2));
					if (magnitude > *xMagnitude || ic == 0) {
						*xMagnitude = magnitude;
						*xOrientation = cvFastArctan(*xYDeriv, *xXDeriv);
					}
				}

				if (*xOrientation >= maxDeg)
					*xOrientation -= maxDeg;
				assert(*xOrientation >= 0 && *xOrientation < 360);
			}
		}

		// in case we have a mask for the image, we have to set the magnitude
		// outside the mask (and at a 1-pixel border within the mask) to zero
		if (imgPyr.getImage(iLevel).hasMask())
			setZeroOutsideMask(gradientMagnitudeTmp, imgPyr.getImage(iLevel).getMask());

		//
		// compute Bin values for each pixel
		//

		// go through all rowsgetN
		for (int iy = 0; iy < levelHeight; ++iy) {
			// initiate pixel pointers
			float* xOrientation = &CV_IMAGE_ELEM(gradientOrientationTmp, float, iy, 0);
			float* xMagnitude = &CV_IMAGE_ELEM(gradientMagnitudeTmp, float, iy, 0);
			std::vector<double*> xHogBins(nBins, NULL);
			for (std::size_t i = 0; i < nBins; ++i)
				xHogBins[i] = &CV_IMAGE_ELEM(hogBins[i].getImage(iLevel), double, iy, 0);

			// go through all pixels in a row
			for (int ix = 0; ix < levelWidth; ++ix, ++xOrientation, ++xMagnitude) {
				// find the correct bin
				float fBin = (*xOrientation) * invBinWidthDeg;
				int iBin = static_cast<int>(roundf(fBin)) % nBins;
				int iBin2 = (fBin - iBin) > 0 ? // the second closest bin for interpolation
						(iBin + 1) % nBins : (iBin - 1 + nBins) % nBins;

				// distribute the magnitude between the two closest bins
				// (use linear interpolation)
				float weight = 1 - std::min<float>(fabs(fBin - iBin), nBins - fBin);
				float weight2 = 1 - weight;
				*(xHogBins[iBin]) = double(weight * (*xMagnitude));
				*(xHogBins[iBin2]) = double(weight2 * (*xMagnitude));

				// increase pointers to hog bin pixels
				for (std::size_t i = 0; i < nBins; ++i) {
					++(xHogBins[i]);
				}
			}
		}

		// compute integral HOG images for each bin
		for (std::size_t i = 0; i < nBins; ++i) {
			assert(hogBins[i].numOfLevels() > iLevel);
			assert(_integralHogBins[i].numOfLevels() > iLevel);
			computeIntegralImage(hogBins[i].getImage(iLevel), _integralHogBins[i].getImage(iLevel));
		}
	}
}

FastHogComputer::VectorType FastHogComputer::getHog(const Box<int>& box, std::size_t iLevel) const
{
	// TODO: method could be rewritten in order to be less expensive
	VectorType hog(_integralHogBins.size());
	fill(hog.begin(), hog.end(), 0);
	iLevel = std::min(iLevel, _integralHogBins[0].numOfLevels() - 1);

	// get the pyramid level we are going to operate on
	int width = _integralHogBins[0].getImage(iLevel)->width;
	int height = _integralHogBins[0].getImage(iLevel)->height;

	// compute the sum of the pixel values for the given box for each bin
	int top = std::min(box.getTop() - 1, height - 1);
	int bottom = std::min(box.getBottom() - 1, height - 1);
	int left = std::min(box.getLeft() - 1, width - 1);
	int right = std::min(box.getRight() - 1, width - 1);
	for (std::size_t i = 0; i < hog.size(); ++i) {
		const IplImageWrapper& img = _integralHogBins[i].getImage(iLevel);
		double sumTopLeft(0), sumTopRight(0), sumBottomLeft(0), sumBottomRight(0);
		if (top >= 0) {
			if (left >= 0)
				sumTopLeft = CV_IMAGE_ELEM(img, double, top, left);
			if (right >= 0)
				sumTopRight = CV_IMAGE_ELEM(img, double, top, right);
		}
		if (bottom >= 0) {
			if (left >= 0)
				sumBottomLeft = CV_IMAGE_ELEM(img, double, bottom, left);
			if (right >= 0)
			sumBottomRight = CV_IMAGE_ELEM(img, double, bottom, right);
		}
		hog[i] = static_cast<ValueType>(sumBottomRight + sumTopLeft - sumBottomLeft - sumTopRight);
	}
	
	// return the hog feature vector
	return hog;
}

FastHogComputer::VectorType FastHogComputer::getIntegralHog(int x, int y, std::size_t iLevel) const {
	VectorType hog(_integralHogBins.size());
	fill(hog.begin(), hog.end(), 0);

	// check whether the point lies inside the image or not
	if (x < 1 || y < 1)
		return hog;
	
	// get the pyramid level we are going to operate on
	iLevel = std::min(iLevel, _integralHogBins[0].numOfLevels() - 1);
	int width = _integralHogBins[0].getImage(iLevel)->width;
	int height = _integralHogBins[0].getImage(iLevel)->height;
	
	// make sure that the point lies in the image
	x = std::min(x - 1, width - 1);
	y = std::min(y - 1, height - 1);
	
	// fill the hog vector
	for (std::size_t i = 0; i < hog.size(); ++i)
		hog[i] = static_cast<ValueType>(CV_IMAGE_ELEM(_integralHogBins[i].getImage(iLevel), double, y, x));
	return hog;
}

std::size_t FastHogComputer::getSizeInBytes() const
{
	std::size_t size(0);
	for (std::size_t i = 0; i < _integralHogBins.size(); ++i)
		for (std::size_t j = 0; j < _integralHogBins[i].numOfLevels(); ++j)
			size += _integralHogBins[i].getImage(j)->height * _integralHogBins[i].getImage(j)->widthStep;
	return size;
}


