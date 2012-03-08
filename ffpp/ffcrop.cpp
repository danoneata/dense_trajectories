#define __STDC_CONSTANT_MACROS 1

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
#include "OutputFormatContext.h"
#include "InputFormat.h"
#include "OutputFormat.h"
#include "Stream.h"
#include "CodecContext.h"
#include "Frame.h"
#include "Packet.h"
#include "ScalerContext.h"
#include "FlipContext.h"

#define DEBUG_FRAMEBLOCK 1000
//#define DEBUG_FRAMEBLOCK 1

// Some utility stuff

using namespace std;
using namespace ffpp;

void usage(char* me)
{
   fprintf(stderr, "Usage: %s [-c <codec>] [-r {dts|pts|sec}] [-s <start>] [-e <end>] [-x <width>] [-y <height>] [-a] [-t <samplerate>] [-C] [-H] [-V] <filename_in> <filename_out>\n"
                   "       -a: resample according to the encoded aspect ratio\n"
                   "       -C: copy frames instead of recoding (conflicts with -c, -x, -y, -a, -H and -V)\n"
                   "       -H: flip each frames horizontally\n"
                   "       -V: flip each frame vertically\n", me);

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


// Here we go

int main(int argc, char *argv[])
{
   // Parse the options

   char* codec = NULL;
   char reference = 'p';
   double frame_start = 0;
   double frame_end = -1;
   int width = 0;
   int height = 0;
   int sample = 1;
   bool aspectFix = false;
   bool copy = false;
   bool horizontalFlip = false;
   bool verticalFlip = false;

   int opt;
   while ((opt = getopt(argc, argv, "c:r:s:e:x:y:t:aCVH")) != -1)
   {
      switch (opt)
      {
         case 'c':
            codec = strdup(optarg);
            break;
         case 'r':
            reference = *optarg;
            break;
         case 's':
            frame_start = atof(optarg);
            break;
         case 'e':
            frame_end = atof(optarg);
            break;
         case 'x':
            width = atoi(optarg);
            break;
         case 'y':
            height = atoi(optarg);
            break;
         case 't':
            sample = atoi(optarg);
            break;
         case 'a':
            aspectFix = true;
            break;
         case 'C':
            copy = true;
            break;
         case 'V':
            verticalFlip = true;
            break;
         case 'H':
            horizontalFlip = true;
            break;
         default:
            usage(argv[0]);
      }
   }

   if (argc < optind+2) usage(argv[0]);
   if (reference != 'd' && reference != 'p' && reference != 's') usage(argv[0]);
   if (copy && (codec || width || height || aspectFix || verticalFlip || horizontalFlip)) usage(argv[0]);

   char* file_in = argv[optind];
   char* file_out = argv[optind+1];

   
   // Init all the crap

   av_register_all();

   
   // Open the input file and get the input format context

   InputFormatContext inFormatCtx(file_in);
   
   dump_format(inFormatCtx.get(), 0, file_in, 0); // Just for debugging

   // Find the input video stream

   Stream inStream = inFormatCtx.GetVideoStream();
   
   // Get the input codec context

   CodecContext inCodecCtx = inStream.GetCodecContext();

   AVRational ratio = inStream.GetAspectRatio();
   if (width == 0 && height == 0) // Same width and height by default
   {
      width = inCodecCtx.GetWidth();
      height = inCodecCtx.GetHeight();
      if (aspectFix)
      {
        if (ratio.num > ratio.den) width = (int)((double)width * ratio.num / ratio.den);
        if (ratio.num < ratio.den) height = (int)((double)height * ratio.den / ratio.num);
      }
   } else { // Try to keep the aspect ratio
      if (width == 0)
      {
        if (aspectFix) width = (int)((double)inCodecCtx.GetWidth() * height / inCodecCtx.GetHeight() * ratio.num / ratio.den);
                  else width = (int)((double)inCodecCtx.GetWidth() * height / inCodecCtx.GetHeight());
      } else if (height == 0) {
        if (aspectFix) height = (int)((double)inCodecCtx.GetHeight() * width / inCodecCtx.GetWidth() * ratio.den / ratio.num);
                  else height = (int)((double)inCodecCtx.GetHeight() * width / inCodecCtx.GetWidth());
      }
   }

   
   // Allocate the output format context, guess the format, set no parameters

   OutputFormat outFormat;
   
   try {
	   outFormat = OutputFormat(0, file_out);
   } catch(std::runtime_error e) {
	   warn("Could not guess the format based on extension, trying the input format");
	   outFormat = OutputFormat(inFormatCtx.GetFormat().GetName(), file_in);
   }
   
   OutputFormatContext outFormatCtx(file_out, outFormat);

   // Create the output video stream

   Stream outStream;
   
   if (codec) outStream = outFormatCtx.AddVideoStream(Codec(codec));
     else outStream = outFormatCtx.AddVideoStream();
   
   outStream.SetFrameRate(inStream.GetFrameRate());
   outStream.SetAspectRatio(inStream.GetAspectRatio());

   // Set the output codec context

   CodecContext outCodecCtx = outStream.GetCodecContext();

   outCodecCtx.SetWidth(width);
   outCodecCtx.SetHeight(height);
   //outCodecCtx.SetBitrate(inCodecCtx.GetBitrate());
   //outCodecCtx.SetBitrate(1000000);

   if (copy)
   {
         outCodecCtx.get()->codec_id = inCodecCtx.get()->codec_id;
         //outCodecCtx.get()->codec_type = inCodecCtx.get()->codec_type;
         outCodecCtx.get()->codec_tag = inCodecCtx.get()->codec_tag;
         //outCodecCtx.get()->has_b_frames = inCodecCtx.get()->has_b_frames;
         //outCodecCtx.get()->bit_rate = inCodecCtx.get()->bit_rate;
         //outCodecCtx.get()->extradata = inCodecCtx.get()->extradata;
         //outCodecCtx.get()->extradata_size = inCodecCtx.get()->extradata_size;
   }

   // Prepare the resampling context
   
   ScalerContext* scalerCtx = 0;

   // Prepare the flip context
   
   FlipContext* flipCtx = 0;
   
   // Open the output file, write the header

   if (outFormatCtx.get()->oformat->flags & AVFMT_NEEDNUMBER)
   {
      if (!av_filename_number_test(outFormatCtx.get()->filename))
         fail("Number token expected in the filename");
   }

   if (!(outFormatCtx.get()->oformat->flags & AVFMT_NOFILE))
   {
      if (url_fopen(&outFormatCtx.get()->pb, outFormatCtx.get()->filename, URL_WRONLY) < 0)
         fail("Could not open the output file");
   }
   
   // This is a convenience hack for correct frame numbering with dts
   if ((reference=='d' || reference=='p') && !strcmp(outFormatCtx.get()->oformat->name, "image2"))
   {
      typedef struct {
          int img_first;
          int img_last;
          int img_number;
          int img_count;
          int is_pipe;
          char path[1024];
      } img2_data_type;

      outFormatCtx.WriteHeader();
      img2_data_type *data = (img2_data_type*)outFormatCtx.get()->priv_data;
      data->img_number = (int)frame_start;
   }

   // Parse the video

   if (outFormatCtx.get()->oformat->flags & AVFMT_RAWPICTURE) fail("Raw dump not supported yet");

   if (frame_start > 0 && reference != 'd')
   {
	  if (reference == 's') inFormatCtx.Seek(frame_start);
	  if (reference == 'p') inFormatCtx.Seek(inStream, (long long)frame_start);
   }

   int haveoutput = 0;
   int havekeyframe = 0;
   long long offset = 0;

   long long dts = 0; // Frame sequence number (Decoding Time Stamp)
   long long pts = 0; // Frame presentation number (Presentation Time Stamp)
   double sec = 0;    // Seconds from the start of the video

   dts = 0;

   boost::optional<Packet> optInPacket;
   while (optInPacket = inFormatCtx.ReadPacket())
   {
	  Packet& inPacket = *optInPacket;
	   
      if (inStream.OwnsPacket(inPacket))
      {
          pts = inPacket.GetPTS();
          sec = inPacket.GetSeconds(inStream);
          
          //fprintf(stderr, "Packet in\n");

          double rts; // Reference timestamp
          if (reference == 'p') rts = pts;
          else if (reference == 's') rts = sec;
          else rts = dts;  	  
          
          boost::optional<Frame> optInFrame;
          if (optInFrame = inCodecCtx.DecodeFrame(inPacket))
          {
        	Frame& inFrame = *optInFrame;
            
            //bool keyframe = inFrame.get()->key_frame;
            bool keyframe = inPacket.get()->flags & PKT_FLAG_KEY;
            
        	const char* ftype = keyframe ? "I" : "P";
	 
            //fprintf(stderr, "Frame in\n");
            
            if (rts >= frame_start && dts%sample == 0)
            {
               if (frame_end >= 0 && rts >= frame_end) break;
	       
               if (!haveoutput)
               {
                  fprintf(stderr, "Starting frame dts=%lld pts=%lld sec=%f (%s)\n", dts, pts, sec, ftype);
                  haveoutput = 1;
                  offset = pts;
               } else if (dts%DEBUG_FRAMEBLOCK == 0 || keyframe) fprintf(stderr, "Copying  frame dts=%lld pts=%lld sec=%f (%s)\n", dts, pts, sec, ftype);
               
               if (keyframe) havekeyframe = 1;
               
               {
                   if (!scalerCtx || !flipCtx)
                   {
                      if (!outCodecCtx.GetEncoder().IsSupportedPixelFormat(inCodecCtx.GetPixelFormat()))
                      {
                         outCodecCtx.SetPixelFormat(outCodecCtx.GetEncoder().GetDefaultPixelFormat());
                      } else {
                         outCodecCtx.SetPixelFormat(inCodecCtx.GetPixelFormat());
                      }

                      dump_format(outFormatCtx.get(), 0, file_out, 1); // Just for debugging

                      // Set the resampling context

                      if (!scalerCtx) scalerCtx = new ScalerContext(inCodecCtx, outCodecCtx);

                      //if (codec) pOutCodec = avcodec_find_encoder_by_name(codec)->id;

                      // Set the flip context

                      if (!flipCtx) flipCtx = new FlipContext(width, height, inCodecCtx.GetPixelFormat(), horizontalFlip, verticalFlip);
                   }

                   Frame outFrame = scalerCtx->Rescale(inFrame);
                   flipCtx->Flip(outFrame);

                   outFrame.SetQuality(inFrame.GetQuality());

                   //fprintf(stderr, "Frame out\n");
               
                   boost::optional<Packet> optOutPacket;
                   if (optOutPacket = outCodecCtx.EncodeFrame(outFrame))
                   {
           	           Packet& outPacket = *optOutPacket;

                       if (copy && havekeyframe)
                       {
                            outStream.OwnPacket(inPacket);
                            inPacket.SetPTS(pts-offset);
                            outFormatCtx.WritePacket(inPacket);
                            //outCodecCtx.get()->frame_number++;

                            //fprintf(stderr, "Packet copy\n");
                       } else {
                           outStream.OwnPacket(outPacket);
                           //outPacket.SetPTS(pts-offset);
                           outFormatCtx.WritePacket(outPacket);
                            
                           //fprintf(stderr, "Packet out\n");
                       }
                   }
               }
            } else {
               if (dts%DEBUG_FRAMEBLOCK == 0) fprintf(stderr, "Skipping frame dts=%lld pts=%lld sec=%f (%s)\n", dts, pts, sec, ftype);
            }

            dts++;
         }
      }
   }
   fprintf(stderr, "Finished frame dts=%lld pts=%lld sec=%f\n", dts-1, pts, sec);

   if (scalerCtx) delete scalerCtx;
   if (flipCtx) delete flipCtx;

   if (codec) free(codec);

   return 0;
}



/*
 * Setting some unusual things, maybe useful someday
 */

//outCodecCtx.get()->gop_size = inCodecCtx.get()->gop_size;
//outCodecCtx.get()->max_b_frames = inCodecCtx.get()->max_b_frames;
//outCodecCtx.get()->flags = inCodecCtx.get()->flags;
//outCodecCtx.get()->max_qdiff = 3;

/*
if (outFormatCtx.get()->oformat->flags & AVFMT_GLOBALHEADER)
   outCodecCtx.get()->flags |= CODEC_FLAG_GLOBAL_HEADER;
*/



/*
 * Workaround for possible missing PTS, seems no longer necessary
 */

//pts_global = packet_in.pts;

/*
uint64_t pts_global = AV_NOPTS_VALUE; // Here we store the current PTS

int pts_get_buffer(struct AVCodecContext *c, AVFrame *pic)
{
   uint64_t *pts = (uint64_t*)av_malloc(sizeof(uint64_t));
   *pts = pts_global;
   pic->opaque = pts; // Here we save the current PTS
   return avcodec_default_get_buffer(c, pic);
}

void pts_release_buffer(struct AVCodecContext *c, AVFrame *pic)
{
   if (pic) av_freep(&pic->opaque);
   avcodec_default_release_buffer(c, pic);
}
*/

//pInCodecCtx->get_buffer = pts_get_buffer;
//pInCodecCtx->release_buffer = pts_release_buffer;

//if (pInFrame->opaque) {
//	pts = offset + *(uint64_t *)pInFrame->opaque;
//}



/*
 * Dunno who needs this for encoding
 */

/*
if (outCodecCtx.get()->coded_frame)
{
   if (outCodecCtx.get()->coded_frame->pts != AV_NOPTS_VALUE)
      outPacket.get()->pts = av_rescale_q(outCodecCtx.get()->coded_frame->pts,
                                    outCodecCtx.get()->time_base, outStream.get()->time_base);
   if (outCodecCtx.get()->coded_frame->key_frame) outPacket.get()->flags |= PKT_FLAG_KEY;
}
*/
