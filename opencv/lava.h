#ifndef opencv_functions_H
#define opencv_functions_H

// LAVA stuff
#include <Images/Image.h>
#include <Images/RGBImage.h>
#include <Images/GrayImage.h>

// OpenCV stuff
#include <opencv/cxcore.h>

lava_ns::RGBImage* iplImg2LavaImgRGB(IplImage *iplImg);
lava_ns::GrayImage* iplImg2LavaGrayImg(IplImage *iplImg);
void iplImg2LavaImgHSV(IplImage *iplImg, lava_ns::GrayImage **hLavaImg, lava_ns::GrayImage **sLavaImg, lava_ns::GrayImage **vLavaImg);
void iplImg2LavaImgYCrCb(IplImage *iplImg, lava_ns::GrayImage **yLavaImg, lava_ns::GrayImage **crLavaImg, lava_ns::GrayImage **cbLavaImg);
template <class T_IMG> T_IMG *loadLavaImg(std::string &filename);
template <> lava_ns::RGBImage *loadLavaImg(std::string &filename);
template <> lava_ns::GrayImage *loadLavaImg(std::string &filename);

#endif
