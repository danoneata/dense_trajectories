# configuration

LAVA2DIR=$(HOME)/c/lava2
LAVA2LIB=/projects/metapkg/opt/lava2/lib

# link .o files with g++
LINK.o=$(CXX) $(LDFLAGS) $(TARGET_ARCH)

# compilation flags
CXXFLAGS_ENV := $(CXXFLAGS)
CXXFLAGS = $(CXXFLAGS_ENV) -I$(HOME)/local/include/ -g -Wall -O2 $(CXXFLAGS_$*)
CXXFLAGS_ffsift = -I$(LAVA2DIR)/include -I$(LAVA2DIR)/src `Magick++-config --cxxflags`

# linking flags
LDFLAGS_ENV := $(LDFLAGS)
LDFLAGS = $(LDFLAGS_ENV) -L$(HOME)/lib $(LDFLAGS_$*)
LDFLAGS_ffsift = -L$(LAVA2LIB) `Magick++-config --ldflags`

# linking libraries
LDLIBS_ENV := $(LDLIBS)
LDLIBS = $(LDLIBS_ENV) -lavutil -lavcodec -lavformat -lswscale $(LDLIBS_$*)
LDLIBS_ffsift = -llava2 `Magick++-config --libs`

# object and program files
OBJS=ScalerContext.o Frame.o Packet.o FormatContext.o InputFormatContext.o OutputFormatContext.o InputFormat.o OutputFormat.o Stream.o CodecContext.o FlipContext.o PixelFormat.o Codec.o
PROG=ffcrop ffcat ffsift

# rules
all: $(PROG)

clean:
	rm -f $(PROG) $(OBJS) ffcrop.o ffsift.o ffcat.o

ffcrop: ffcrop.o $(OBJS) 

ffcat: ffcat.o $(OBJS) 

ffsift: ffsift.o $(OBJS)

