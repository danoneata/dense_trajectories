// standard stuff
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <cassert>
#include <string>

// STL stuff
//#include <vector>
//#include <list>


// LAVA stuff
#include <Images/Image.h>
#include <Images/RGBImage.h>
#include <Images/GrayImage.h>

// OpenCV stuff
#include <opencv/cv.h>
#include <opencv/highgui.h>


/**
 * converts an IplImage to a Lava GrayImage
 */
lava_ns::GrayImage *iplImg2LavaGrayImg(IplImage *iplImg)
{
	using namespace lava_ns;

	int width = iplImg->width;
	int height = iplImg->height;

	// convert color to gray and gray to floating point
	IplImage *imgGray8U = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, 1);
	IplImage *imgGray = cvCreateImage(cvSize(width, height), IPL_DEPTH_64F, 1);
	cvCvtColor(iplImg, imgGray8U, CV_RGB2GRAY);
	cvConvertScale(imgGray8U, imgGray);

	// prepare a LAVA version of the IplImage
	lava_ns::GrayImage *imgLava = new lava_ns::GrayImage(width, height);

	// copy row for row
	for (int iRow = 0; iRow < height; iRow++)
		memcpy((*imgLava)[iRow], imgGray->imageData + imgGray->widthStep * iRow, width * imgGray->depth / 8);

// 	for (int x = 0; x < imgGray8U->width; x++)
// 		for (int y = 0; y < imgGray8U->height; y++)
// 			CV_IMAGE_ELEM(imgGray8U, uchar, y, x) = (uchar)imgLava->getpixel(x, y);
// 	cvNamedWindow("testImage", CV_WINDOW_AUTOSIZE);
// 	cvShowImage("testImage", imgGray8U);
// 	cvWaitKey(0);
// 	cvDestroyWindow("testImage");

	// clean up
	cvReleaseImage(&imgGray8U);
	cvReleaseImage(&imgGray);

	return imgLava;
}


void iplImg2LavaImgHSV(IplImage *iplImg, lava_ns::GrayImage **hLavaImg, lava_ns::GrayImage **sLavaImg, lava_ns::GrayImage **vLavaImg)
{
	using namespace lava_ns;

	int width = iplImg->width;
	int height = iplImg->height;

	// convert color to gray and gray to floating point
	IplImage *imgHsv8U = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, 3);
	IplImage *imgHsv = cvCreateImage(cvSize(width, height), IPL_DEPTH_64F, 3);
	cvCvtColor(iplImg, imgHsv8U, CV_RGB2HSV);
	cvConvert(imgHsv8U, imgHsv);

	// prepare a LAVA version of the IplImage
	*hLavaImg = new lava_ns::GrayImage(width, height);
	*sLavaImg = new lava_ns::GrayImage(width, height);
	*vLavaImg = new lava_ns::GrayImage(width, height);

	// copy position for position
//	double *pH, *pS, *pV, *pHSV;
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			(*hLavaImg)->setpixel(x, y, CV_IMAGE_ELEM(imgHsv, double, y, x * imgHsv->nChannels) * 255.0 / 180.0);
			(*sLavaImg)->setpixel(x, y, CV_IMAGE_ELEM(imgHsv, double, y, x * imgHsv->nChannels + 1));
			(*vLavaImg)->setpixel(x, y, CV_IMAGE_ELEM(imgHsv, double, y, x * imgHsv->nChannels + 2));
		}
	}

//  	IplImage *imgGray8U = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, 1);
// 	uchar tmp;
// 	cvNamedWindow("testImage1", CV_WINDOW_AUTOSIZE);
// 	cvNamedWindow("testImage2", CV_WINDOW_AUTOSIZE);
// 	cvNamedWindow("testImage3", CV_WINDOW_AUTOSIZE);
// 	for (int x = 0; x < imgGray8U->width; x++)
// 		for (int y = 0; y < imgGray8U->height; y++)
// 			CV_IMAGE_ELEM(imgGray8U, uchar, y, x) = (uchar)(*hLavaImg)->getpixel(x, y);
// 	cvShowImage("testImage1", imgGray8U);
// 	for (int x = 0; x < imgGray8U->width; x++)
// 		for (int y = 0; y < imgGray8U->height; y++)
// 			CV_IMAGE_ELEM(imgGray8U, uchar, y, x) = (uchar)(*sLavaImg)->getpixel(x, y);
// 	cvShowImage("testImage2", imgGray8U);
// 	for (int x = 0; x < imgGray8U->width; x++)
// 		for (int y = 0; y < imgGray8U->height; y++)
// 			CV_IMAGE_ELEM(imgGray8U, uchar, y, x) = (uchar)(*vLavaImg)->getpixel(x, y);
// 	cvShowImage("testImage3", imgGray8U);
// 	cvWaitKey(0);
// 	cvDestroyWindow("testImage1");
// 	cvDestroyWindow("testImage2");
// 	cvDestroyWindow("testImage3");
// 	cvReleaseImage(&imgGray8U);

	// clean up
	cvReleaseImage(&imgHsv);
	cvReleaseImage(&imgHsv8U);
}

lava_ns::RGBImage *iplImg2LavaImgRGB(IplImage *iplImg)
{
	using namespace lava_ns;

	int width = iplImg->width;
	int height = iplImg->height;

	// prepare a LAVA version of the IplImage
	lava_ns::RGBImage *lavaImg = new lava_ns::RGBImage(width, height);

	// copy position for position
	RGBPixel rgb;
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			rgb.r = (unsigned short)MAX(MIN(round(CV_IMAGE_ELEM(iplImg, char, y, x * iplImg->nChannels + 2)), 255), 0);
			rgb.g = (unsigned short)MAX(MIN(round(CV_IMAGE_ELEM(iplImg, char, y, x * iplImg->nChannels + 1)), 255), 0);
			rgb.b = (unsigned short)MAX(MIN(round(CV_IMAGE_ELEM(iplImg, char, y, x * iplImg->nChannels + 0)), 255), 0);
			lavaImg->setpixel(x, y, rgb);
		}
	}

	return lavaImg;
}

void iplImg2LavaImgYCrCb(IplImage *iplImg, lava_ns::GrayImage **yLavaImg, lava_ns::GrayImage **crLavaImg, lava_ns::GrayImage **cbLavaImg)
{
	using namespace lava_ns;
	using namespace std;

	int width = iplImg->width;
	int height = iplImg->height;

	// convert color to gray and gray to floating point
	IplImage *imgYCrCb8U = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, 3);
	IplImage *imgYCrCb = cvCreateImage(cvSize(width, height), IPL_DEPTH_32F, 3);
	IplImage *imgY = cvCreateImage(cvSize(width, height), IPL_DEPTH_32F, 1);
	IplImage *imgYHighPassed = cvCreateImage(cvSize(width, height), IPL_DEPTH_32F, 1);
	cvCvtColor(iplImg, imgYCrCb8U, CV_BGR2YCrCb);
	cvConvert(imgYCrCb8U, imgYCrCb);
	cvSplit(imgYCrCb, imgY, NULL, NULL, NULL);
	cvSmooth(imgY, imgYHighPassed, CV_GAUSSIAN, 2 * ((width / 5) / 2) + 1, 2 * ((height / 5) / 2) + 1);
// 	cvSmooth(imgY, imgYHighPassed, CV_GAUSSIAN, 15);
// 	for (int y = 0; y < height; y++)
// 		for (int x = 0; x < width; x++)
// 			cout << CV_IMAGE_ELEM(imgY, float, y, x) << " - " << CV_IMAGE_ELEM(imgYHighPassed, float, y, x) << endl;
	cvSub(imgY, imgYHighPassed, imgYHighPassed);

	// make sure that the value range is between 0 and 255
	double maxVal, minVal;
	cvMinMaxLoc(imgYHighPassed, &minVal, &maxVal);
	cvConvertScale(imgYHighPassed, imgYHighPassed, 255.0f / (maxVal - minVal), -255.0 * minVal / (maxVal - minVal));

	// prepare a LAVA version of the IplImage
	*yLavaImg = new lava_ns::GrayImage(width, height);
	*crLavaImg = new lava_ns::GrayImage(width, height);
	*cbLavaImg = new lava_ns::GrayImage(width, height);

	// copy position for position
	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++) {
// cout << CV_IMAGE_ELEM(imgYHighPassed, float, y, x) << endl;
			(*yLavaImg)->setpixel(x, y, CV_IMAGE_ELEM(imgYHighPassed, float, y, x));
			(*crLavaImg)->setpixel(x, y, CV_IMAGE_ELEM(imgYCrCb, float, y, x * imgYCrCb->nChannels + 1));
			(*cbLavaImg)->setpixel(x, y, CV_IMAGE_ELEM(imgYCrCb, float, y, x * imgYCrCb->nChannels + 2));
		}

// 	IplImage *imgGray8U = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, 1);
// 	cvNamedWindow("testImage1", CV_WINDOW_AUTOSIZE);
// 	cvNamedWindow("testImage2", CV_WINDOW_AUTOSIZE);
// 	cvNamedWindow("testImage3", CV_WINDOW_AUTOSIZE);
// 	for (int x = 0; x < imgGray8U->width; x++)
// 		for (int y = 0; y < imgGray8U->height; y++)
// 			CV_IMAGE_ELEM(imgGray8U, uchar, y, x) = (uchar)(*yLavaImg)->getpixel(x, y);
// 	cvShowImage("testImage1", imgGray8U);
// 	for (int x = 0; x < imgGray8U->width; x++)
// 		for (int y = 0; y < imgGray8U->height; y++)
// 			CV_IMAGE_ELEM(imgGray8U, uchar, y, x) = (uchar)(*crLavaImg)->getpixel(x, y);
// 	cvShowImage("testImage2", imgGray8U);
// 	for (int x = 0; x < imgGray8U->width; x++)
// 		for (int y = 0; y < imgGray8U->height; y++)
// 			CV_IMAGE_ELEM(imgGray8U, uchar, y, x) = (uchar)(*cbLavaImg)->getpixel(x, y);
// 	cvShowImage("testImage3", imgGray8U);
// 	cvWaitKey(0);
// 	cvDestroyWindow("testImage1");
// 	cvDestroyWindow("testImage2");
// 	cvDestroyWindow("testImage3");
// 	cvReleaseImage(&imgGray8U);

	// clean up
	cvReleaseImage(&imgYCrCb);
	cvReleaseImage(&imgY);
	cvReleaseImage(&imgYHighPassed);
	cvReleaseImage(&imgYCrCb8U);
}

template <class T_IMG>
T_IMG *loadLavaImg(std::string &filename)
{
	return NULL;
}

template <>
lava_ns::GrayImage *loadLavaImg(std::string &filename)
{
	// precondition
	assert(filename.length() > 0);

	IplImage *loadedImg = cvLoadImage(filename.c_str());
	lava_ns::GrayImage *lavaImg;
	if (NULL == loadedImg)
		lavaImg = new lava_ns::GrayImage();
	else {
		lavaImg = iplImg2LavaGrayImg(loadedImg);
		if (NULL == lavaImg)
			lavaImg = new lava_ns::GrayImage();
		cvReleaseImage(&loadedImg);
	}

	// postcondition
	assert(NULL != lavaImg);
	return lavaImg;
}

template <>
lava_ns::RGBImage *loadLavaImg(std::string &filename)
{
	// precondition
	assert(filename.length() > 0);

	IplImage *loadedImg = cvLoadImage(filename.c_str());
	lava_ns::RGBImage *lavaImg;
	if (NULL == loadedImg)
		lavaImg = new lava_ns::RGBImage();
	else {
		lavaImg = iplImg2LavaImgRGB(loadedImg);
		if (NULL == lavaImg)
			lavaImg = new lava_ns::RGBImage();
		cvReleaseImage(&loadedImg);
	}

	// postcondition
	assert(NULL != lavaImg);
	return lavaImg;
}
