#define __STDC_CONSTANT_MACROS 1
#define __LAVA2_NO_ITS 1

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

extern "C" {
	#include <libavcodec/avcodec.h>
	#include <libavformat/avformat.h>
}

#include <stdexcept>

#include "InputFormatContext.h"
#include "InputFormat.h"
#include "Stream.h"
#include "CodecContext.h"
#include "Frame.h"
#include "Packet.h"
#include "ScalerContext.h"

#include <Descriptors/CornerDescriptor.h>
#include <Images/ImagePyramid.h>
#include <Images/GrayImage.h>
#include <Detectors/HarrisLaplace.h>
#include <Descriptors/SIFT.h>

#define DEBUG_FRAMEBLOCK 25
//#define DEBUG_FRAMEBLOCK 1

// Some utility stuff

using namespace std;
using namespace ffpp;

void usage(char* me)
{
   fprintf(stderr, "Usage: %s [-s <start>{dts|pts|sec}] [-e <end>{dts|pts|sec}] [-t <every>{dts|pts|sec}] [-x <width>] [-y <height>] <filename_in>\n", me);

   exit(-1);
}

void fail(char* message)
{
   fprintf(stderr, "Error: %s\n", message);

   exit(-2);
}

void warn(char* message)
{
   fprintf(stderr, "Warning: %s\n", message);
}

void debug(char* message)
{
   fprintf(stdout, "%s\n", message);
}

typedef lava_ns::GrayImage Image;

class FramePointReceiver: public lava_ns::PointReceiver
{
	protected:
		int width;
		int height;
		double timebase;
		long long ts;
		lava_ns::Descriptor<Image>& descriptor;
	
	public:
		FramePointReceiver(int width, int height, double timebase, lava_ns::Descriptor<Image>& descriptor):
			width(width), height(height), timebase(timebase), ts(0), descriptor(descriptor)
		{
			std::cout << "x y t x-norm y-norm t-sec scale sift(128)" << std::endl;
		}
		
		void setts(long long ts)
		{
			this->ts = ts;
		}
		
		void insert(lava_ns::CornerDescriptor & cd)
		{
			descriptor.describe(cd);
			
			int i;

			std::cout << cd.getx()              << " " << cd.gety()               << " " << ts          << " ";
			std::cout << (float)cd.getx()/width << " " << (float)cd.gety()/height << " " << ts*timebase << " ";
			std::cout << cd.getscale();
			for (i=0;i<cd.size();++i) std::cout << " " << cd[i];
			std::cout << std::endl;
		}
};

// Here we go

struct Timestamp
{
	Timestamp(): set(false)
	{}
	
	Timestamp(const char* str): set(false)
	{
		char* end;
		
		value = strtod(str, &end);
		if (end != str)
		{
			set = true;
			type = *end;
		}
	}
	
	bool set;
	double value;
	char type;
};

int main(int argc, char *argv[])
{
   // Parse the options

   Timestamp frame_start;
   Timestamp frame_end;
   Timestamp frame_every;
   int width = 0;
   int height = 0;

   int opt;
   while ((opt = getopt(argc, argv, "s:e:t:x:y:")) != -1)
   {
      switch (opt)
      {
         case 's':
            frame_start = Timestamp(optarg);
            break;
         case 'e':
            frame_end = Timestamp(optarg);
            break;
         case 't':
            frame_every = Timestamp(optarg);
            break;
         case 'x':
            width = atoi(optarg);
            break;
         case 'y':
            height = atoi(optarg);
            break;
         default:
            usage(argv[0]);
      }
   }

   if (argc < optind+1) usage(argv[0]);
   if (  frame_end.set &&   frame_end.type != 'd' &&   frame_end.type != 'p' &&   frame_end.type != 's') usage(argv[0]);
   if (frame_start.set && frame_start.type != 'd' && frame_start.type != 'p' && frame_start.type != 's') usage(argv[0]);
   if (frame_every.set && frame_every.type != 'd' && frame_every.type != 'p' && frame_every.type != 's') usage(argv[0]);

   char* file_in = argv[optind];

   // Init all the crap

   av_register_all();

   
   // Open the input file and get the input format context

   InputFormatContext inFormatCtx(file_in);
   
   dump_format(inFormatCtx.get(), 0, file_in, 0); // Just for debugging

   // Find the input video stream

   Stream inStream = inFormatCtx.GetVideoStream();
   
   // Get the input codec context

   CodecContext inCodecCtx = inStream.GetCodecContext();

   if (width == 0 && height == 0) // Same width and height by default
   {
      width = inCodecCtx.GetWidth();
      height = inCodecCtx.GetHeight();
   } else { // Try to keep the aspect ratio
      if (width == 0) width = (int)((float)inCodecCtx.GetWidth() * height/inCodecCtx.GetHeight());
      else if (height == 0) height = (int)((float)inCodecCtx.GetHeight() * width/inCodecCtx.GetWidth());
   }
   
   // Set the resampling context
   
   ScalerContext* pScalerCtx = 0;
   if (inCodecCtx.GetPixelFormat() >= 0)
	pScalerCtx = new ScalerContext(inCodecCtx, width, height, PIX_FMT_GRAY8);

   // Set detector/descriptor stuff

   Image img(width, height);
   lava_ns::ImagePyramid<Image> pyramid;
   pyramid.setScaleFactor(1.2);
   lava_ns::SIFT<Image> descriptor;
   descriptor.setImagePyramid(&pyramid);
   descriptor.setMainParam(61,12,4,8);
   FramePointReceiver receiver(width, height, av_q2d(inStream.GetTimeBase()), descriptor);
   lava_ns::HarrisLaplace<Image> detector;
   detector.setPosThreshold(300);
   detector.setReceiver(&receiver);
   detector.setImage(&pyramid);
   
   // Parse the video

   if (frame_start.set)
   {
	  if (frame_start.type == 's') inFormatCtx.Seek(frame_start.value);
	  if (frame_start.type == 'p') inFormatCtx.Seek(inStream, frame_start.value);	  
   }

   int haveoutput = 0;

   long long dts = 0; // Frame sequence number (Decoding Time Stamp)
   long long pts = 0; // Frame presentation number (Presentation Time Stamp)
   double sec = 0;    // Seconds from the start of the video

   boost::optional<Packet> optInPacket;
   while (optInPacket = inFormatCtx.ReadPacket())
   {
	  Packet& inPacket = *optInPacket;
	   
      if (inStream.OwnsPacket(inPacket))
      {
          pts = inPacket.GetPTS();
          sec = inPacket.GetSeconds(inStream);

          boost::optional<Frame> optInFrame;
          if (optInFrame = inCodecCtx.DecodeFrame(inPacket))
          {
        	Frame& inFrame = *optInFrame;
        	  
        	const char* ftype = inFrame.get()->key_frame ? "I" : "P";
	 
        	bool read = true;
        	if (frame_start.set)
        	{
        		read = false;
        		if (frame_start.type == 'd' && dts >= frame_start.value) read = true;
        		if (frame_start.type == 'p' && pts >= frame_start.value) read = true;
        		if (frame_start.type == 's' && sec >= frame_start.value) read = true;	
        		if (read) frame_start.set = false;
        	}
        	
    		if (read && frame_every.set)
    		{
   				if (frame_start.type != frame_every.type)
   				{
    				if (frame_every.type == 'd') frame_start.value = dts;
    				if (frame_every.type == 'p') frame_start.value = pts;
    				if (frame_every.type == 's') frame_start.value = sec;
   					frame_start.type = frame_every.type; 
   				}
   				frame_start.set = true;
   				frame_start.value += frame_every.value; 
    		}
        	
            if (read)
            {
               if (frame_end.set)
               {
            	   if (frame_end.type == 'd' && dts >= frame_end.value) break;
            	   if (frame_end.type == 'p' && pts >= frame_end.value) break;
            	   if (frame_end.type == 's' && sec >= frame_end.value) break;
               }
	       
               if (!haveoutput)
               {
                  fprintf(stderr, "Starting frame dts=%lld pts=%lld sec=%f (%s)\n", dts, pts, sec, ftype);
                  haveoutput = 1;
               } //else if (dts%DEBUG_FRAMEBLOCK == 0 || inFrame.get()->key_frame) fprintf(stderr, "Copying  frame dts=%lld pts=%lld sec=%f (%s)\n", dts, pts, sec, ftype);

	       if (!pScalerCtx) pScalerCtx = new ScalerContext(inCodecCtx, width, height, PIX_FMT_GRAY8);
               Frame outFrame = pScalerCtx->Rescale(inFrame);
            
               int x,y;
		       for (y=0;y<height;++y)
		       {
		    	   for (x=0;x<width;++x)
		    	   {
		    		   img[y][x] = outFrame.Get8Pixel(x,y);
	               }
		       }

		       receiver.setts(pts);
               pyramid.setImage(&img);
               detector.detect();
            } else {
               //if (dts%DEBUG_FRAMEBLOCK == 0) fprintf(stderr, "Skipping frame dts=%lld pts=%lld sec=%f (%s)\n", dts, pts, sec, ftype);
            }

            dts++;
         }
      }
   }
   fprintf(stderr, "Finished frame dts=%lld pts=%lld sec=%f\n", dts-1, pts, sec);

   if (pScalerCtx) delete pScalerCtx;

   return 0;
}
