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

// Some utility stuff

using namespace std;
using namespace ffpp;

void usage(char* me)
{
   fprintf(stderr, "Usage: %s <filename_in> [...] <filename_out>\n", me);
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

   int opt;
   while ((opt = getopt(argc, argv, "")) != -1)
   {
      switch (opt)
      {
         default:
            usage(argv[0]);
      }
   }

   if (argc < optind+2) usage(argv[0]);

   char* file_in = argv[optind];
   char* file_out = argv[argc-1];
   
   // Init all the crap

   av_register_all();
   
   // Open the input file and get the input format context

   InputFormatContext inFormatCtx(file_in);
   
   dump_format(inFormatCtx.get(), 0, file_in, 0); // Just for debugging

   // Find the input video stream

   Stream inStream = inFormatCtx.GetVideoStream();
   
   // Get the input codec context

   CodecContext inCodecCtx = inStream.GetCodecContext();

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

   Stream outStream = outFormatCtx.AddVideoStream();
   
   outStream.SetFrameRate(inStream.GetFrameRate());
   outStream.SetAspectRatio(inStream.GetAspectRatio());

   // Set the output codec context

   CodecContext outCodecCtx = outStream.GetCodecContext();

   outCodecCtx.SetWidth(inCodecCtx.GetWidth());
   outCodecCtx.SetHeight(inCodecCtx.GetHeight());
   //outCodecCtx.SetBitrate(inCodecCtx.GetBitrate());

   outCodecCtx.get()->codec_id = inCodecCtx.get()->codec_id;
   //outCodecCtx.get()->codec_type = inCodecCtx.get()->codec_type;
   outCodecCtx.get()->codec_tag = inCodecCtx.get()->codec_tag;
   //outCodecCtx.get()->has_b_frames = inCodecCtx.get()->has_b_frames;
   //outCodecCtx.get()->bit_rate = inCodecCtx.get()->bit_rate;
   //outCodecCtx.get()->extradata = inCodecCtx.get()->extradata;
   //outCodecCtx.get()->extradata_size = inCodecCtx.get()->extradata_size;

   //if (codec) pOutCodec = avcodec_find_encoder_by_name(codec)->id;
   
   // Open the output file, write the header

   if (!(outFormatCtx.get()->oformat->flags & AVFMT_NOFILE))
   {
      if (url_fopen(&outFormatCtx.get()->pb, outFormatCtx.get()->filename, URL_WRONLY) < 0)
         fail("Could not open the output file");
   }
   
   // Parse the video

   if (outFormatCtx.get()->oformat->flags & AVFMT_RAWPICTURE) fail("Raw dump not supported yet");

   long long globaloffset = 0;

   while (optind < argc-1)
   {
        char* file_in = argv[optind];
        fprintf(stderr, "Appending%s\n", file_in);
        
        InputFormatContext inFormatCtx(file_in);
        Stream inStream = inFormatCtx.GetVideoStream();
        
        long long localoffset = -1;
        long long pts = 0; // Frame presentation number (Presentation Time Stamp)

        boost::optional<Packet> optInPacket;
        while (optInPacket = inFormatCtx.ReadPacket())
        {
            Packet& inPacket = *optInPacket;
	   
            if (inStream.OwnsPacket(inPacket))
            {
                pts = inPacket.GetPTS();
                
                if (localoffset < 0) localoffset = pts;
            
                //fprintf(stderr, "dts=%lu, pts=%lu\n", inPacket.get()->dts, inPacket.get()->pts);

                outStream.OwnPacket(inPacket);
                inPacket.SetPTS(pts-localoffset+globaloffset);
                outFormatCtx.WritePacket(inPacket);
            }
        }
        
        globaloffset += pts-localoffset + 1;
        optind++;
   }

   return 0;
}
