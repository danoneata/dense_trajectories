#ifndef FEATTRACK_H_
#define FEATTRACK_H_

// OpenCV stuff
#include <opencv/cv.h>
#include <opencv/highgui.h>

#include <opencv/Video.h>
#include <opencv/IplImageWrapper.h>
#include <opencv/IplImagePyramid.h>
#include <opencv/FastHogComputer.h>

#include <geometry/Box.hpp>
#include <geometry/Size.hpp>

#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <list>
#include <cmath>

#include <boost/numeric/ublas/io.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

using std::string;
using std::vector;

bool fullOrientation = true;
float epsilon = 0.05;
float scaleFactor = sqrt(2);
const float min_flow = 0.4*0.4; // TBS 
const float PI = 3.14159;

// parameters for tracking
int win_size = 10;
double quality = 0.001;
double min_distance = 5;
const float min_var = sqrt(3);
const float max_var = 50;
const float max_dis = 20;

// parameters for multiscale
int scale_flag = 1; // 0: singlescale; 1: multiscale TBS
int scale_num = 8;  // For singlescale, set scale_num to 1u TBS
int scale_show = 0; // show tracks at which scale TBS
const float scale_stride = sqrt(2);

int normalize_tag = 0;  // 0: normalize L2 norm to 1; 1: normalize length to 1 TBS
const int max_frameNum = 1000000;
const float mindist_merge = 2.0;
const float mindist_remove = 3.0;
// parameters for shot boundary
int shot_tag = 0; // 0: keep boundary tracks; 1: remove boundary tracks
const int shot_margin = 1;  // TBS

typedef struct Sequence
{
  int width;
  int height;
  int length;
}Sequence;

typedef struct TrackerInfo
{
  int trackLength; /*length of the trajectory*/
  int initGap; /*initial gap for feature detection*/
  int norm; /*0: length normalization; 1: L1; 2: L2*/
}TrackerInfo;

typedef struct DescInfo
{
  int nBins; /*number of bins for vector quantization*/
  int fullOrientation; /*0: 180 degree; 1: 360 degree*/
  int norm; /*1: L1 normalization; 2: L2 normalization*/
  float threshold; /*threshold for normalization*/
  int flagThre; /* whether thresholding or not */
  int nxCells; /*number of cells in x direction*/
  int nyCells;
  int ntCells;
  int dim;
  int blockHeight; /*size of the block for computing the descriptor*/
  int blockWidth;
  int flag; /*0: don't compute and output the descriptor; 1: otherwise*/
}DescInfo; 

typedef struct DescMat
{
  int height;
  int width;
  int nBins;
  float* desc;
}DescMat;

class PointDesc
{
  public:
    std::vector<float> hog;
    std::vector<float> hof;
    std::vector<float> mbhX;
    std::vector<float> mbhY;
    CvPoint2D32f point;

    PointDesc(const DescInfo& hogInfo, const DescInfo& hofInfo, const DescInfo& mbhInfo, const CvPoint2D32f& point_)
      : hog(hogInfo.nxCells * hogInfo.nyCells * hogInfo.nBins),
      hof(hofInfo.nxCells * hofInfo.nyCells * hofInfo.nBins),
      mbhX(mbhInfo.nxCells * mbhInfo.nyCells * mbhInfo.nBins),
      mbhY(mbhInfo.nxCells * mbhInfo.nyCells * mbhInfo.nBins),
      point(point_)
  { }
};

class Track
{
  public:
    std::list<PointDesc> pointDescs;
    int maxNPoints;

    Track(int maxNPoints_)
      : maxNPoints(maxNPoints_)
    {}

    void addPointDesc(const PointDesc& point)
    {
      pointDescs.push_back(point);
      if (pointDescs.size() > maxNPoints + 2) {
        pointDescs.pop_front();
        //printf("the track is too long!\n");
      }
    }
};
#endif /*FEATTRACK_H_*/
