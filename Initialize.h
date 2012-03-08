#ifndef INITIALIZE_H_
#define INITIALIZE_H_

#include "FeatTrack.h"

void InitTrackerInfo(TrackerInfo* tracker, int track_length, int init_gap)
{
	tracker->trackLength = track_length;
	tracker->initGap = init_gap;
    tracker->norm = 0;
}

DescMat* InitDescMat(int height, int width, int nBins)
{
    DescMat* descMat = (DescMat*)malloc(sizeof(DescMat));
	descMat->height = height;
	descMat->width = width;
	descMat->nBins = nBins;
	descMat->desc = (float*)malloc(height*width*nBins*sizeof(float));
	memset( descMat->desc, 0, height*width*nBins*sizeof(float));
    return descMat;
}

void ReleDescMat( DescMat* descMat)
{
	free(descMat->desc);
	free(descMat);
}

void InitDescInfo(DescInfo* descInfo, int nBins, int flag, int orientation, int size, int nxy_cell, int nt_cell)
{
	descInfo->nBins = nBins;
	descInfo->fullOrientation = orientation;
	descInfo->norm = 2;
	descInfo->threshold = min_flow;
	descInfo->flagThre = flag;
	descInfo->nxCells = nxy_cell;
	descInfo->nyCells = nxy_cell;
	descInfo->ntCells = nt_cell;
	descInfo->dim = descInfo->nBins*descInfo->nxCells*descInfo->nyCells;
	descInfo->blockHeight = size;
	descInfo->blockWidth = size;
	descInfo->flag = 1;
}

void InitSequence(char* video, Sequence& seq)
{
	Video capture(video);
	// get the number of frames in the video
	int counter = 0;
	while( true ) {
		IplImageWrapper frame = 0;
        frame = capture.getFrame();
		
        if( !frame )
            break;
		counter++;
	  
		if (!capture.nextFrame())
			break;
    }
	//std::cerr << "the video length: " << counter << std::endl;
	seq.length = counter;
	seq.width = capture.getWidth();
	seq.height = capture.getHeight();
	// relocate the video to the first frame
	long long frameIndex = 0;
	capture.seek(frameIndex);
}

#endif /*INITIALIZE_H_*/
