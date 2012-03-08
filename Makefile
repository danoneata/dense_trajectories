# opencv2.2 needed + setmyenv lear no_epd no_locvideo 2.2
#
# set the binaries that have to be built
TARGETS := FeatTrack

# set the build configuration set 
BUILD := release
#BUILD := debug
BIT := 64

# set bin and build dirs
BUILDDIR := .build_$(BUILD)$(BIT)
BINDIR := $(BUILD)$(BIT)

ALEXDIR := /home/clear/gaidon/progs/alex
LOCDIR := $(HOME)/local
# opencv < 2.3 required (or problem with ffpp)
OPENCV_VER := 2.2

# include directories
INCLUDEDIRS := \
 	$(LOCDIR)/opencv/opencv$(OPENCV_VER)/include \
	/usr/include/ffmpeg
#	$(LOCDIR)/video/include \

# library directories
LIBDIRS := \
 	$(LOCDIR)/opencv/opencv$(OPENCV_VER)/lib
#	$(LOCDIR)/video/lib \

# libraries (without the prefix "lib") 
LIBS := \
	boost_program_options boost_regex boost_system boost_filesystem \
	opencv_core opencv_highgui opencv_video opencv_imgproc \
	avformat avdevice avutil avcodec swscale
#	lapack cblas atlas blas \

# set which libraries are used by which executable
LDLIBS = $(addprefix -L, $(LIBDIRS)) $(addprefix -l, $(LIBS))

# set some flags and compiler/linker specific commands
CXXFLAGS_debug := -ggdb -D_SHOW_TRACKS
CXXFLAGS_release := -O3 -DNDEBUG -ggdb -DN_SHOW_TRACKS
CXXFLAGS = -m$(BIT) -pipe -D __STDC_CONSTANT_MACROS -D STD=std -Wall $(CXXFLAGS_$(BUILD)) -I. -I/opt/include $(addprefix -I, $(INCLUDEDIRS))

LDFLAGS_debug := -ggdb
LDFLAGS_release := -O3 -ggdb
LDFLAGS = -m$(BIT) -L/opt/lib -pipe -Wall $(LDFLAGS_$(BUILD)) $(addprefix -L, $(LIBDIRS))

include make/generic.mk
