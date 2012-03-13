#include "FeatTrack.h"
#include "Descriptors.h"
#include "Initialize.h"

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <vector>

//std namespaces
//namespace std { using namespace __gnu_cxx; }
using std::string;
using std::vector;
using namespace cv;

IplImageWrapper image, prev_image, grey, prev_grey;
IplImagePyramid grey_pyramid, prev_grey_pyramid, eig_pyramid;

float* fscales = 0; /*float scale values*/

// structure used to save computed descriptors to disk in a binary
// format (.siftgeo) that is readable by bigimbaz:
//   1st the geom structure,
//   then the dimension of the descriptor (int)
//   then the descriptor (dim floats)
typedef struct {
  struct {
    // geom header as in bigimbaz's siftgeo format
    float x, y, t, xvar, yvar, length, fscale, dummy1, dummy2;
  } geom;
  int dim;
  float * descriptor;
} point_t;

int usage() {
  fprintf(stderr,"usage:\n\
      FeatTrack <video_file> <output_file>.siftgeo (if '0', print descriptor instead of saving it)\n\
                [<track_length_in_frames> <stride_in_pixels> (default: 15 5)] \n\
                [<beginning_frame> <end_frame> (default: begin and end of video)]\n\
                [<descriptor name: hog, hof, mbh, track[r] or all[r]> (default: all)]\n\
                [<refresh_rate (default: 1)>]\n");
  return -1;
}

// Functions for splitting string by a separator.
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    return split(s, delim, elems);
}

// Checks if the given frame xx is within one
// of the given ranges (low[ii], high[ii]).
int is_in_range(int xx, vector<int> low, vector<int> high) {
  int len = low.size();
  for (int ii = 0; ii < len; ii++) {
    if ((low[ii] <= xx) && (xx <= high[ii])) {
      return 1;
    }
  }
  return 0;
}

int main( int argc, char** argv )
{
  int return_raw_track = 0;
  int start_frame = 1;
  int end_frame = max_frameNum;
  int frameNum = 0;
  TrackerInfo tracker;
  DescInfo hogInfo;
  DescInfo hofInfo;
  DescInfo mbhInfo;
  string desc_name = "all";
  int init_gap = 1;
  int track_length = 15;
  float min_distance = 5.;
  if( argc < 3 )
    return usage();
  int patch_size = 32;
  int nxy_cell = 2;	
  int nt_cell = 3;
  Video capture(argv[1]);
  char* out_filename = argv[2];	
  int SW_STDOUT = 0;
  FILE * fo;
  vector<int> start_frames;
  vector<int> end_frames;

  if (strcmp(out_filename, "0") == 0)
    // Don't save the descriptor; pass it as STDOUT.
    SW_STDOUT = 1;

  Sequence sequence;
  InitSequence(argv[1], sequence);

  /*read in parameters*/
  if( argc > 4 ) {
    track_length = atoi(argv[3]);	
    min_distance = atof(argv[4]);
  }
  if( argc > 6 ) {
    // start_frame = atoi(argv[5]);
    // end_frame = atoi(argv[6]);
    std::vector<std::string> str_start_frames = split(argv[5], '_');
    std::vector<std::string> str_end_frames = split(argv[6], '_');
    int nr_start_frames = str_start_frames.size();
    int nr_end_frames = str_end_frames.size();
    assert(nr_start_frames == nr_end_frames);
    int nr_frames = nr_start_frames;
    for(int ii = 0; ii < nr_frames; ii++) {
      start_frames.push_back(atoi(str_start_frames[ii].c_str()));
      end_frames.push_back(atoi(str_end_frames[ii].c_str()));
    }
  }
  else {
    // If no arguments specified for start_frame and end_frame
    // then use defaults.
    start_frames.push_back(start_frame);
    end_frames.push_back(end_frame);
  }
  if( argc > 7 ) {
    desc_name = argv[7];
    if (desc_name == "allr") {
        // use raw trajectory descriptor
        desc_name = "all";
        return_raw_track = 1;
    } else if (desc_name == "trackr") {
        // use raw trajectory descriptor
        desc_name = "track";
        return_raw_track = 1;
    } else if (desc_name != "all" && desc_name != "hog" && desc_name != "hof" && desc_name != "mbh" && desc_name != "track")
        return usage();
  }
  if( argc > 8 ) {
    init_gap = atoi(argv[8]);
    if( init_gap <= 0 )
      return usage();
  }

  InitTrackerInfo(&tracker, track_length, init_gap);
  InitDescInfo(&hogInfo, 8, 0, 1, patch_size, nxy_cell, nt_cell);
  InitDescInfo(&hofInfo, 9, 1, 1, patch_size, nxy_cell, nt_cell);
  InitDescInfo(&mbhInfo, 8, 0, 1, patch_size, nxy_cell, nt_cell);

  start_frame -= 1; /*for KTH, frameNum starts with 1*/
  end_frame -= 1;

  if (!SW_STDOUT) {
    // open the output file
    fo = fopen (out_filename, "w");
    assert (fo);
  }

  // total number of points tracked
  int nbpts = 0;

  // the descriptor to be saved
  point_t desc;
  desc.geom.dummy1 = 0.0; // dummy field to have correct geom size
  desc.geom.dummy2 = 0.0; // dummy field to have correct geom size

  // compute the dimensionality of the full vector
  desc.dim = 0;
  if(desc_name == "all" || desc_name == "track")
    desc.dim += 2 * tracker.trackLength;           // track descriptor (x and y)
  if(desc_name == "all" || desc_name == "hog")
    desc.dim += hogInfo.ntCells * hogInfo.dim;     // hog descriptor
  if(desc_name == "all" || desc_name == "hof")
    desc.dim += hofInfo.ntCells * hofInfo.dim;     // hof descriptor
  if(desc_name == "all" || desc_name == "mbh")
    desc.dim += 2 * mbhInfo.ntCells * mbhInfo.dim; // mbhx + mbhy descriptor

#ifdef _SHOW_TRACKS
  cvNamedWindow( "LkDemo", 0 );
#endif

  std::vector<std::list<Track> > xyScaleTracks;
  int init_counter = 0; /*indicate when to initialize the tracker*/
  while( true ) {
    IplImageWrapper frame = 0;
    unsigned int i;

    frame = capture.getFrame();
    frameNum = capture.getFrameIndex();

    if( !frame )
      break;
    // if( frameNum >= start_frame && frameNum <= end_frame ) {
    if (is_in_range(frameNum, start_frames, end_frames)) {
      if( !image ) {
        /*allocate all the buffers*/
        image = IplImageWrapper( cvGetSize(frame), 8, 3 );
        image->origin = frame->origin;
        prev_image= IplImageWrapper( cvGetSize(frame), 8, 3 );
        prev_image->origin = frame->origin;
        grey = IplImageWrapper( cvGetSize(frame), 8, 1 );
        grey_pyramid = IplImagePyramid( cvGetSize(frame), 8, 1, scale_stride );
        prev_grey = IplImageWrapper( cvGetSize(frame), 8, 1 );
        prev_grey_pyramid = IplImagePyramid( cvGetSize(frame), 8, 1, scale_stride );
        eig_pyramid = IplImagePyramid( cvGetSize(frame), 32, 1, scale_stride );

        cvCopy( frame, image, 0 );
        cvCvtColor( image, grey, CV_BGR2GRAY );
        grey_pyramid.rebuild(grey);

        scale_num = std::min<std::size_t>(scale_num, grey_pyramid.numOfLevels());
        fscales = (float*)cvAlloc(scale_num*sizeof(float)); 
        xyScaleTracks.resize(scale_num);

        for( int ixyScale = 0; ixyScale < scale_num; ++ixyScale ) {
          std::list<Track>& tracks = xyScaleTracks[ixyScale];
          fscales[ixyScale] = pow(scale_stride, ixyScale);

          /*find good features at each scale*/
          IplImage *grey_temp = 0, *eig_temp = 0;
          std::size_t temp_level = (std::size_t)ixyScale;
          grey_temp = cvCloneImage(grey_pyramid.getImage(temp_level));
          eig_temp = cvCloneImage(eig_pyramid.getImage(temp_level));
          std::vector<CvPoint2D32f> points(0);	
          cvDenseSample(grey_temp, eig_temp, points, quality, min_distance);
          /*save the feature points*/
          for( i = 0; i < points.size(); i++ ) {
            Track track(tracker.trackLength);
            PointDesc point(hogInfo, hofInfo, mbhInfo, points[i]);
            track.addPointDesc(point);
            tracks.push_back(track);
          }

          cvReleaseImage( &grey_temp );
          cvReleaseImage( &eig_temp );
        }
      }

      cvCopy( frame, image, 0 );
      cvCvtColor( image, grey, CV_BGR2GRAY );
      grey_pyramid.rebuild(grey);

      // printf("frameNum: %d\n", frameNum);
      // std::cerr << "frameNum: " << frameNum << std::endl;
      if( frameNum > 0 ) {
        init_counter++;
        for( int ixyScale = 0; ixyScale < scale_num; ++ixyScale ) {
          std::vector<CvPoint2D32f> points_in(0);
          std::list<Track>& tracks = xyScaleTracks[ixyScale];
          for (std::list<Track>::iterator iTrack = tracks.begin(); iTrack != tracks.end(); ++iTrack) {
            CvPoint2D32f point = iTrack->pointDescs.back().point;
            points_in.push_back(point); 
          }
          int count = points_in.size();
          /*track points_in at each scale*/
          IplImage *prev_grey_temp = 0, *grey_temp = 0;
          std::size_t temp_level = ixyScale;
          prev_grey_temp = cvCloneImage(prev_grey_pyramid.getImage(temp_level));
          grey_temp = cvCloneImage(grey_pyramid.getImage(temp_level));

          Mat prev_grey_mat = cvarrToMat(prev_grey_temp);
          Mat grey_mat = cvarrToMat(grey_temp);

          std::vector<int> status(count);
          std::vector<CvPoint2D32f> points_out(count);

          // compute the dense optical flow field
          IplImage* flow = cvCreateImage(cvGetSize(grey_temp), IPL_DEPTH_32F, 2);
          Mat flow_mat = cvarrToMat(flow);
          calcOpticalFlowFarneback( prev_grey_mat, grey_mat, flow_mat,
              sqrt(2)/2.0, 5, 10, 2, 7, 1.5, OPTFLOW_FARNEBACK_GAUSSIAN );

          // track the points by optical flow interpolation
          OpticalFlowTracker(flow, points_in, points_out, status);			

          // compute the different integral images
          int width = grey_temp->width;
          int height = grey_temp->height;
          DescMat* hogMat = InitDescMat(height, width, hogInfo.nBins);
          HogComp(prev_grey_temp, hogMat, hogInfo);

          DescMat* hofMat = InitDescMat(height, width, hofInfo.nBins);
          HofComp(flow, hofMat, hofInfo);

          DescMat* mbhMatX = InitDescMat(height, width, mbhInfo.nBins);
          DescMat* mbhMatY = InitDescMat(height, width, mbhInfo.nBins);
          MbhComp(flow, mbhMatX, mbhMatY, mbhInfo);	

          // compute the descriptors from the above integral images
          i = 0;
          for (std::list<Track>::iterator iTrack = tracks.begin(); iTrack != tracks.end(); ++i) {
            if( status[i] == 1 ) {
              PointDesc& pointDesc = iTrack->pointDescs.back(); 	
              CvPoint2D32f prev_point = points_in[i];
              CvScalar rect;
              getRect(&rect, prev_point, cvSize(width, height), hogInfo);
              Box<double> box_temp( rect.val[0], rect.val[1], rect.val[2], rect.val[3] );
              FastHogComputer::VectorType hog_temp(hogInfo.dim);
              hog_temp = getDesc(hogMat, box_temp, hogInfo);
              for( int m  = 0; m < hogInfo.dim; m++ )
                pointDesc.hog[m] = hog_temp[m];

              FastHogComputer::VectorType hof_temp(hofInfo.dim);
              hof_temp = getDesc(hofMat, box_temp, hofInfo);
              for( int m  = 0; m < hofInfo.dim; m++ )
                pointDesc.hof[m] = hof_temp[m];

              FastHogComputer::VectorType mbh_temp(mbhInfo.dim);
              mbh_temp = getDesc(mbhMatX, box_temp, mbhInfo);
              for( int m  = 0; m < mbhInfo.dim; m++ )
                pointDesc.mbhX[m] = mbh_temp[m];

              mbh_temp = getDesc(mbhMatY, box_temp, mbhInfo);
              for( int m  = 0; m < mbhInfo.dim; m++ )
                pointDesc.mbhY[m] = mbh_temp[m];
              PointDesc point(hogInfo, hofInfo, mbhInfo, points_out[i]);
              iTrack->addPointDesc(point);
#ifdef _SHOW_TRACKS
              // draw this track
              std::list<PointDesc>& descs = iTrack->pointDescs;
              std::list<PointDesc>::iterator iDesc = descs.begin();
              float ll = descs.size();
              CvPoint2D32f point0 = iDesc->point; 
              float jj = 0;
              point0.x *= fscales[ixyScale]; // map the point to first scale
              point0.y *= fscales[ixyScale];

              for (iDesc++; iDesc != descs.end(); ++iDesc, ++jj) {
                CvPoint2D32f point1 = iDesc->point;
                point1.x *= fscales[ixyScale]; 
                point1.y *= fscales[ixyScale];

                if (fabs(point1.x - point0.x) + fabs(point1.y - point0.y) > 0.5*min_distance) {
                  // don't redraw
                  cvLine(image, cvPointFrom32f(point0),
                    //cvPointFrom32f(point1), CV_RGB(0,0,255), 2, 8,0);
                    cvPointFrom32f(point1), CV_RGB(0,cvFloor(255.0*(jj+1.0)/ll),0), 2, 8,0);
                }
                point0 = point1;
              }
              //cvCircle(image, cvPointFrom32f(point0), 2, CV_RGB(255,0,0), -1, 8,0);
#endif
              ++iTrack;
            }
            else 
              iTrack = tracks.erase(iTrack);
          }
          ReleDescMat(hogMat);
          ReleDescMat(hofMat);
          ReleDescMat(mbhMatX);
          ReleDescMat(mbhMatY);
          cvReleaseImage( &prev_grey_temp );
          cvReleaseImage( &grey_temp );
          cvReleaseImage( &flow );
        }	    	

        for( int ixyScale = 0; ixyScale < scale_num; ++ixyScale ) {
          std::list<Track>& tracks = xyScaleTracks[ixyScale];
          for (std::list<Track>::iterator iTrack = tracks.begin(); iTrack != tracks.end(); ) {
            if( iTrack->pointDescs.size() >= tracker.trackLength+1 ) {
              std::vector<CvPoint2D32f> trajectory(tracker.trackLength+1);
              std::list<PointDesc>& descs = iTrack->pointDescs;
              std::list<PointDesc>::iterator iDesc = descs.begin();

              for (int count = 0; count <= tracker.trackLength; ++iDesc, ++count) {
                trajectory[count].x = iDesc->point.x*fscales[ixyScale];
                trajectory[count].y = iDesc->point.y*fscales[ixyScale];
              }
              float mean_x(0), mean_y(0), var_x(0), var_y(0), length(0);
              // check if track is valid
              int is_valid;
              if (return_raw_track == 1) {
                // make a copy to leave trajectory unchanged to return the real trajectory
                std::vector<CvPoint2D32f> copy_trajectory(trajectory);
                is_valid = isValid(copy_trajectory, mean_x, mean_y, var_x, var_y, length);
                if (is_valid == 1) {
                  // remove last point (like in isValid) to get proper length
                  trajectory.pop_back();
                  // printf("%d %d %d\n", (int) tracker.trackLength, (int) trajectory.size(), (int) copy_trajectory.size());
                }
              }
              else {
                // don't  make a copy: modify trajectory in-place to contain
                // it's shape descriptor (normalized displacement vectors)
                is_valid = isValid(trajectory, mean_x, mean_y, var_x, var_y, length);
              }
              if( is_valid == 1 ) {
                nbpts++;
                // print the track information
                //printf("%d\t", frameNum);
                if (SW_STDOUT)
                  printf("%f\t%f\t%d\t", mean_x, mean_y, frameNum);
                //printf("%f\t%f\t", var_x, var_y);
                //printf("%f\t", length);
                //printf("%f\t", fscales[ixyScale]);
                desc.geom.x = mean_x;
                desc.geom.y = mean_y;
                desc.geom.t = frameNum;
                desc.geom.xvar = var_x;
                desc.geom.yvar = var_y;
                desc.geom.length = length;
                desc.geom.fscale = fscales[ixyScale];

                // allocate space for the descriptor
                desc.descriptor = (float *) malloc (desc.dim*sizeof(float));

                // current index in the full descriptor
                int idesc = 0;

                // print the track descriptor
                if(desc_name == "all" || desc_name == "track")
                  for (int count = 0; count < tracker.trackLength; ++count)	{
                    if (SW_STDOUT)
                      printf("%f\t%f\t", trajectory[count].x,trajectory[count].y );
                    desc.descriptor[idesc] = trajectory[count].x;
                    idesc++;
                    desc.descriptor[idesc] = trajectory[count].y;
                    idesc++;
                  }

                int t_stride;
                // print the hog descriptor
                if(desc_name == "all" || desc_name == "hog") {
                  iDesc = descs.begin();
                  t_stride = cvFloor(tracker.trackLength/hogInfo.ntCells);
                  for( int n = 0; n < hogInfo.ntCells; n++ ) {
                    float vec[hogInfo.dim];
                    for( int m = 0; m < hogInfo.dim; m++ )
                      vec[m] = 0;

                    for( int t = 0; t < t_stride; t++, iDesc++ ) {
                      for( int m = 0; m < hogInfo.dim; m++ )
                        vec[m] += iDesc->hog[m];
                    }
                    for( int m = 0; m < hogInfo.dim; m++ ) {
                      if (SW_STDOUT)
                        printf("%f\t", vec[m]/float(t_stride));
                      desc.descriptor[idesc] = vec[m]/float(t_stride);
                      idesc++;
                    }
                  }
                }

                // print the hof descriptor
                if(desc_name == "all" || desc_name == "hof") {
                  iDesc = descs.begin();
                  t_stride = cvFloor(tracker.trackLength/hofInfo.ntCells);
                  for( int n = 0; n < hofInfo.ntCells; n++ ) {
                    float vec[hofInfo.dim];
                    for( int m = 0; m < hofInfo.dim; m++ )
                      vec[m] = 0;

                    for( int t = 0; t < t_stride; t++, iDesc++ ) {
                      for( int m = 0; m < hofInfo.dim; m++ )
                        vec[m] += iDesc->hof[m];
                    }
                    for( int m = 0; m < hofInfo.dim; m++ ) {
                      if (SW_STDOUT)
                        printf("%f\t", vec[m]/float(t_stride));
                      desc.descriptor[idesc] = vec[m]/float(t_stride);
                      idesc++;
                    }
                  }
                }

                if(desc_name == "all" || desc_name == "mbh") {
                  // print the MBH(x) descriptor
                  iDesc = descs.begin();
                  t_stride = cvFloor(tracker.trackLength/mbhInfo.ntCells);
                  for( int n = 0; n < mbhInfo.ntCells; n++ ) {
                    float vec[mbhInfo.dim];
                    for( int m = 0; m < mbhInfo.dim; m++ )
                      vec[m] = 0;

                    for( int t = 0; t < t_stride; t++, iDesc++ ) {
                      for( int m = 0; m < mbhInfo.dim; m++ )
                        vec[m] += iDesc->mbhX[m];
                    }
                    for( int m = 0; m < mbhInfo.dim; m++ ) {
                      if (SW_STDOUT)
                        printf("%f\t", vec[m]/float(t_stride));
                      desc.descriptor[idesc] = vec[m]/float(t_stride);
                      idesc++;
                    }
                  }

                  // print the MBH(y) descriptor
                  iDesc = descs.begin();
                  t_stride = cvFloor(tracker.trackLength/mbhInfo.ntCells);
                  for( int n = 0; n < mbhInfo.ntCells; n++ ) {
                    float vec[mbhInfo.dim];
                    for( int m = 0; m < mbhInfo.dim; m++ )
                      vec[m] = 0;

                    for( int t = 0; t < t_stride; t++, iDesc++ ) {
                      for( int m = 0; m < mbhInfo.dim; m++ )
                        vec[m] += iDesc->mbhY[m];
                    }
                    for( int m = 0; m < mbhInfo.dim; m++ ) {
                      if (SW_STDOUT)
                        printf("%f\t", vec[m]/float(t_stride));
                      desc.descriptor[idesc] = vec[m]/float(t_stride);
                      idesc++;
                    }
                  }
                }

                // write the descriptor to file
                if (SW_STDOUT) {
                  printf("\n");
                  fflush(stdout);
                }
                else {
                  fwrite (&desc.geom , sizeof (desc.geom), 1, fo);
                  fwrite (&desc.dim, sizeof (desc.dim), 1, fo);
                  fwrite (desc.descriptor, sizeof(float), desc.dim, fo);
                }
                free(desc.descriptor);

              }
              
              iTrack = tracks.erase(iTrack);
            }
            else
              iTrack++;
          }
        }	

        if( init_counter == tracker.initGap ) { /*initialize every initGap frames*/
          init_counter = 0;
          for (int ixyScale = 0; ixyScale < scale_num; ++ixyScale) {
            std::list<Track>& tracks = xyScaleTracks[ixyScale];			
            std::vector<CvPoint2D32f> points_in(0);
            std::vector<CvPoint2D32f> points_out(0);
            for(std::list<Track>::iterator iTrack = tracks.begin(); iTrack != tracks.end(); iTrack++, i++) {
              std::list<PointDesc>& descs = iTrack->pointDescs; /*the last point in the quenue*/
              CvPoint2D32f point = descs.back().point;
              points_in.push_back(point);
            }

            IplImage *grey_temp = 0, *eig_temp = 0;
            std::size_t temp_level = (std::size_t)ixyScale;
            grey_temp = cvCloneImage(grey_pyramid.getImage(temp_level));
            eig_temp = cvCloneImage(eig_pyramid.getImage(temp_level));

            cvDenseSample(grey_temp, eig_temp, points_in, points_out, quality, min_distance);

            for( i = 0; i < points_out.size(); i++) {
              Track track(tracker.trackLength);
              PointDesc point(hogInfo, hofInfo, mbhInfo, points_out[i]);
              track.addPointDesc(point);
              tracks.push_back(track);
            }
            cvReleaseImage( &grey_temp );
            cvReleaseImage( &eig_temp );
          }
        }
      }

      cvCopy( frame, prev_image, 0 );
      cvCvtColor( prev_image, prev_grey, CV_BGR2GRAY );
      prev_grey_pyramid.rebuild(prev_grey);
    }

#ifdef _SHOW_TRACKS
    cvShowImage( "LkDemo", image);
    unsigned int c;
    c = cvWaitKey(50);   
    if( (char)c == 27 )
      break;
#endif

    if (!capture.nextFrame())
      break;
  }

  if (!SW_STDOUT) {
    // close the output file
    fclose (fo);
    fprintf (stderr, "Found %d tracks of dimension %d\n", nbpts, desc.dim);
  }

#ifdef _SHOW_TRACKS
  cvDestroyWindow("LkDemo");
#endif


  return 0;
}

#ifdef _EiC
main(1,"lkdemo.c");
#endif
