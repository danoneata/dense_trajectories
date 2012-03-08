#ifndef DESCRIPTORS_H_
#define DESCRIPTORS_H_

#include "FeatTrack.h"

/* the rectangle for computing the descriptor */
void getRect(CvScalar* rect, // returned rectangle
    CvPoint2D32f point, // interest point position
    CvSize size, // the size of the image
    DescInfo descInfo) // parameters about the descriptor
{
  int x_min = descInfo.blockWidth/2;
  int y_min = descInfo.blockHeight/2;
  int x_max = size.width - descInfo.blockWidth;
  int y_max = size.height - descInfo.blockHeight;

  CvPoint2D32f point_temp;

  float temp = point.x - x_min;
  temp = (temp > 0) ? temp : 0;
  temp = (temp > x_max) ? x_max : temp;
  point_temp.x = temp;

  temp = point.y - y_min;
  temp = (temp > 0) ? temp : 0;
  temp = (temp > y_max) ? y_max : temp;
  point_temp.y = temp;

  rect->val[0] = point_temp.x;
  rect->val[1] = point_temp.y;
  rect->val[2] = descInfo.blockWidth;
  rect->val[3] = descInfo.blockHeight;
}

/* compute integral histogram for the whole image */
void BuildDescMat(IplImage* xComp, // x gradient or flow component
    IplImage* yComp, // y gradient or flow component
    DescMat* descMat, // output integral histogram
    DescInfo descInfo) // parameters about the descriptor
{
  float fullAngle = descInfo.fullOrientation ? 360 : 180;
  int nBins = descInfo.flagThre ? descInfo.nBins-1 : descInfo.nBins;
  float angleBase = fullAngle/float(nBins);
  int width = descMat->width;
  int height = descMat->height;
  int histDim = descMat->nBins;
  float *sum = (float*)cvAlloc(histDim*sizeof(float));
  int index = 0;
  for(int i = 0; i < height; i++) {
    const float* xcomp = (const float*)(xComp->imageData + xComp->widthStep*i);
    const float* ycomp = (const float*)(yComp->imageData + yComp->widthStep*i);
    for(int j = 0; j < histDim; j++)
      sum[j] = 0;
    for(int j = 0; j < width; j++, index++) {
      float shiftX = xcomp[j];
      float shiftY = ycomp[j];
      float magnitude0 = sqrt(shiftX*shiftX+shiftY*shiftY);
      float magnitude1 = magnitude0;
      int bin0, bin1;
      if(descInfo.flagThre == 1 && magnitude0 <= descInfo.threshold) {
        bin0 = nBins;
        magnitude0 = 1.0;
        bin1 = 0;
        magnitude1 = 0;
      }
      else {
        float orientation = cvFastArctan(shiftY, shiftX);

        if(orientation > fullAngle)
          orientation -= fullAngle;
        //orientation = descInfo.fullOrientation ? orientation : orientation-180.;

        float fbin = orientation/angleBase;
        bin0 = static_cast<int>(roundf(fbin)) % nBins;
        bin1 = (fbin - bin0) > 0 ? (bin0 + 1) % nBins 
          : (bin0 - 1 + nBins) % nBins;

        float weight0 = 1 - std::min<float>(fabs(fbin - bin0), nBins - fbin);
        float weight1 = 1 - weight0;

        magnitude0 *= weight0;
        magnitude1 *= weight1;
      }

      sum[bin0] += magnitude0;
      sum[bin1] += magnitude1;

      int temp0 = index*descMat->nBins;
      if(i == 0) {
        for(int m = 0; m < descMat->nBins; m++)
          descMat->desc[temp0++] = sum[m];
      }
      else {
        int temp1 = (index - width)*descMat->nBins;
        for(int m = 0; m < descMat->nBins; m++)
          descMat->desc[temp0++] = descMat->desc[temp1++]+sum[m];
      }
    }
  }
  cvFree(&sum);
}

/* extract a descriptor from the integral histogram */
FastHogComputer::VectorType getDesc(DescMat* descMat, // input integral histogram
    Box<double> box, // rectangle area for the descriptor
    DescInfo descInfo) // parameters about the descriptor
{
  int descDim = descInfo.dim;
  int height = descMat->height;
  int width = descMat->width;

  FastHogComputer::VectorType vec(descDim);
  int xStride = box.getWidth()/descInfo.nxCells;
  int yStride = box.getHeight()/descInfo.nyCells;
  int xOffset = box.getX();
  int yOffset = box.getY();
  CvPoint pt;

  /* iterate over different cells */
  int iDesc = 0;
  for (int iX = 0; iX < descInfo.nxCells; ++iX)
    for (int iY = 0; iY < descInfo.nyCells; ++iY) {
      int left = xOffset + iX*xStride - 1;  
      int right = left + xStride;
      right = (right < width-1) ? right : width-1;
      int top = yOffset + iY*yStride - 1;
      int bottom = top + yStride;
      bottom = (bottom < height-1) ? bottom : height-1;

      int TopLeft = top*width+left;
      int TopRight = top*width+right;
      int BottomLeft = bottom*width+left;
      int BottomRight = bottom*width+right;
      TopLeft *= descInfo.nBins;
      TopRight *= descInfo.nBins;
      BottomLeft *= descInfo.nBins;
      BottomRight *= descInfo.nBins;

      for (int i = 0; i < descInfo.nBins; ++i, ++iDesc) {
        double sumTopLeft(0), sumTopRight(0), sumBottomLeft(0), sumBottomRight(0);
        if (top >= 0) {
          if (left >= 0)
            sumTopLeft = descMat->desc[TopLeft+i];
          if (right >= 0)
            sumTopRight = descMat->desc[TopRight+i];
        }
        if (bottom >= 0) {
          if (left >= 0)
            sumBottomLeft = descMat->desc[BottomLeft+i];
          if (right >= 0)
            sumBottomRight = descMat->desc[BottomRight+i];
        }
        float temp = sumBottomRight + sumTopLeft
          - sumBottomLeft - sumTopRight;
        /* why temp can be negative? */
        temp = (temp > 0) ? temp : 0;
        temp += epsilon;
        vec[iDesc] = temp;
      }
    }
  if (descInfo.norm == 1) /* L1 normalization */
    vec *= 1 / boost::numeric::ublas::norm_1(vec);
  else /* L2 normalization */
    vec *= 1 / boost::numeric::ublas::norm_2(vec);

  return vec;
}

void HogComp(IplImage* img, DescMat* descMat, DescInfo descInfo)
{
  int width = descMat->width;
  int height = descMat->height;
  IplImage* imgX = cvCreateImage(cvSize(width,height), IPL_DEPTH_32F, 1);
  IplImage* imgY = cvCreateImage(cvSize(width,height), IPL_DEPTH_32F, 1);
  cvSobel(img, imgX, 1, 0, 1);
  cvSobel(img, imgY, 0, 1, 1);
  BuildDescMat(imgX, imgY, descMat, descInfo);
  cvReleaseImage(&imgX);
  cvReleaseImage(&imgY);
}

void HofComp(IplImage* flow, DescMat* descMat, DescInfo descInfo)
{
  int width = descMat->width;
  int height = descMat->height;
  IplImage* xComp = cvCreateImage(cvSize(width, height), IPL_DEPTH_32F, 1);
  IplImage* yComp = cvCreateImage(cvSize(width, height), IPL_DEPTH_32F, 1);
  for(int i = 0; i < height; i++) {
    const float* f = (const float*)(flow->imageData + flow->widthStep*i);
    float* xf = (float*)(xComp->imageData + xComp->widthStep*i);
    float* yf = (float*)(yComp->imageData + yComp->widthStep*i);
    for(int j = 0; j < width; j++) {
      xf[j] = f[2*j];
      yf[j] = f[2*j+1];
    }
  }
  BuildDescMat(xComp, yComp, descMat, descInfo);
  cvReleaseImage(&xComp);
  cvReleaseImage(&yComp);
}

void MbhComp(IplImage* flow, DescMat* descMatX, DescMat* descMatY, DescInfo descInfo)
{
  int width = descMatX->width;
  int height = descMatX->height;
  IplImage* flowX = cvCreateImage(cvSize(width,height), IPL_DEPTH_32F, 1);
  IplImage* flowY = cvCreateImage(cvSize(width,height), IPL_DEPTH_32F, 1);
  IplImage* flowXdX = cvCreateImage(cvSize(width,height), IPL_DEPTH_32F, 1);
  IplImage* flowXdY = cvCreateImage(cvSize(width,height), IPL_DEPTH_32F, 1);
  IplImage* flowYdX = cvCreateImage(cvSize(width,height), IPL_DEPTH_32F, 1);
  IplImage* flowYdY = cvCreateImage(cvSize(width,height), IPL_DEPTH_32F, 1);

  /* extract the x and y components of the flow */
  for(int i = 0; i < height; i++) {
    const float* f = (const float*)(flow->imageData + flow->widthStep*i);
    float* fX = (float*)(flowX->imageData + flowX->widthStep*i);
    float* fY = (float*)(flowY->imageData + flowY->widthStep*i);
    for(int j = 0; j < width; j++) {
      fX[j] = 100*f[2*j];
      fY[j] = 100*f[2*j+1];
      //printf("x:%f y:%f\t", f[2*j], f[2*j+1]);
    }
    //printf("\n");
  }

  cvSobel(flowX, flowXdX, 1, 0, 1);
  cvSobel(flowX, flowXdY, 0, 1, 1);
  cvSobel(flowY, flowYdX, 1, 0, 1);
  cvSobel(flowY, flowYdY, 0, 1, 1);

  BuildDescMat(flowXdX, flowXdY, descMatX, descInfo);
  BuildDescMat(flowYdX, flowYdY, descMatY, descInfo);

  cvReleaseImage(&flowX);	
  cvReleaseImage(&flowY);
  cvReleaseImage(&flowXdX);
  cvReleaseImage(&flowXdY);
  cvReleaseImage(&flowYdX);
  cvReleaseImage(&flowYdY);
}

/* tracking interesting points by interpolating in the optical field */
void OpticalFlowTracker(IplImage* flow, // the optical field
    std::vector<CvPoint2D32f>& points_in, // input interesting point positions
    std::vector<CvPoint2D32f>& points_out, // output interesting point positions
    std::vector<int>& status) // flag for successfully tracked or not
{
  if(points_in.size() != points_out.size())
    std::cerr << "points numbers don't match!" << std::endl;
  if(points_in.size() != status.size())
    std::cerr << "the status number doesn't match the point number" << std::endl;
  int width = flow->width;
  int height = flow->height;
  const int num = points_in.size();

  for(int i = 0; i < num; i++) {
    CvPoint2D32f point_in = points_in[i];
    std::list<float> xs;
    std::list<float> ys; 
    int x = cvFloor(point_in.x);
    int y = cvFloor(point_in.y);	
    for(int m = x-1; m <= x+1; m++)
      for(int n = y-1; n <= y+1; n++) {
        int p = m;
        p = (p > 0) ? p : 0;
        p = (p < width-1) ? p : width-1;
        int q = n;
        q = (q > 0) ? q : 0;
        q = (q < height-1) ? q : height-1;
        const float* f = (const float*)(flow->imageData + flow->widthStep*q);
        xs.push_back(f[2*p]);
        ys.push_back(f[2*p+1]);
      }

    xs.sort();
    ys.sort();
    int size = xs.size();
    size /= 2;
    for(int m = 0; m < size; m++) {
      xs.pop_back();
      ys.pop_back();
    }

    CvPoint2D32f offset;
    offset.x = xs.back();
    offset.y = ys.back();
    xs.clear();
    ys.clear();
    CvPoint2D32f point_out;
    point_out.x = point_in.x + offset.x;
    point_out.y = point_in.y + offset.y;
    points_out[i] = point_out;
    //printf("the %dth point: x:%f\ty:%f\n", i, point_out.x, point_out.y);
    if( point_out.x > 0 && point_out.x < width && point_out.y > 0 && point_out.y < height )
      status[i] = 1;
    else
      status[i] = -1;		 
  }
}	

/* check whether a trajectory is valid or not */
int isValid(std::vector<CvPoint2D32f>& track, float& mean_x, float& mean_y, float& var_x, float& var_y, float& length) 
{
  int size = track.size();
  for(int i = 0; i < size; i++) {
    mean_x += track[i].x;
    mean_y += track[i].y;
  }
  mean_x /= size;
  mean_y /= size;

  for(int i = 0; i < size; i++) {
    track[i].x -= mean_x;
    var_x += track[i].x*track[i].x;
    track[i].y -= mean_y;
    var_y += track[i].y*track[i].y;
  }
  var_x /= size;
  var_y /= size;
  var_x = sqrt(var_x);	
  var_y = sqrt(var_y);	

  if(var_x < min_var && var_y < min_var)
    return 0;

  if( var_x > max_var || var_y > max_var ) 
    return 0;

  for(int i = 1; i < size; i++) {
    float temp_x = track[i].x - track[i-1].x;
    float temp_y = track[i].y - track[i-1].y;
    length += sqrt(temp_x*temp_x+temp_y*temp_y);
    track[i-1].x = temp_x;
    track[i-1].y = temp_y;
  }

  float len_thre = length*0.7;
  for( int i = 0; i < size-1; i++ ) {
    float temp_x = track[i].x;
    float temp_y = track[i].y;
    float temp_dis = sqrt(temp_x*temp_x + temp_y*temp_y);
    if( temp_dis > max_dis && temp_dis > len_thre ) 
      return 0;
  }

  track.pop_back();
  for(int i = 0; i < size-1; i++) {
    track[i].x /= length;
    track[i].y /= length;
  }
  return 1;
}

void cvDenseSample(IplImage* grey, IplImage* eig, std::vector<CvPoint2D32f>& points,
    const double quality, const double min_distance)
{
  int width = cvFloor(grey->width/min_distance);
  int height = cvFloor(grey->height/min_distance);
  double maxVal = 0;
  cvCornerMinEigenVal(grey, eig, 3, 3);
  cvMinMaxLoc(eig, 0, &maxVal, 0, 0, 0);
  const double threshold = maxVal*quality;

  int offset = cvFloor(min_distance/2);
  for(int i = 0; i < height; i++) {
    for(int j = 0; j < width; j++) {
      int x = cvFloor(j*min_distance+offset);
      int y = cvFloor(i*min_distance+offset);
      if(CV_IMAGE_ELEM(eig, float, y, x) > threshold) {
        CvPoint2D32f point;
        point.x = x;
        point.y = y;
        points.push_back(point);
      }
    }
  }
}

void cvDenseSample(IplImage* grey, IplImage* eig, std::vector<CvPoint2D32f>& points_in, 
    std::vector<CvPoint2D32f>& points_out, const double quality, const double min_distance)
{
  int width = cvFloor(grey->width/min_distance);
  int height = cvFloor(grey->height/min_distance);
  int area = width*height;
  double maxVal = 0;
  cvCornerMinEigenVal(grey, eig, 3, 3);
  cvMinMaxLoc(eig, 0, &maxVal, 0, 0, 0);
  const double threshold = maxVal*quality;

  std::vector<int> counters(area);
  for(int i = 0; i < points_in.size(); i++) {
    CvPoint2D32f point = points_in[i];
    if(point.x >= min_distance*width || point.y >= min_distance*height)
      continue;
    int x = cvFloor(point.x/min_distance);
    int y = cvFloor(point.y/min_distance);
    counters[y*width+x]++;
  }

  int index = 0;
  int offset = cvFloor(min_distance/2);
  for(int i = 0; i < height; i++) {
    for(int j = 0; j < width; j++, index++) {
      if(counters[index] == 0) {
        int x = cvFloor(j*min_distance+offset);
        int y = cvFloor(i*min_distance+offset);
        if(CV_IMAGE_ELEM(eig, float, y, x) > threshold) {
          CvPoint2D32f point;
          point.x = x;
          point.y = y;
          points_out.push_back(point);
        }
      }
    }
  }
}

#endif /*DESCRIPTORS_H_*/
