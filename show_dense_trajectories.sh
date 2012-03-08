#!/usr/bin/env bash

source $HOME/.profile

video="$1"
shift

L="$1"        # max length of trajectory in frames
shift

W="$1"        # stride in pixels for dense sampling
shift

R="$1"        # refresh rate in frames to resample points
shift

# default values
if [ "$L" = "" ]; then
  L=15
fi
if [ "$W" = "" ]; then
  W=5
fi
if [ "$R" = "" ]; then
  W=1
fi

desc="mbh"

exe="/home/clear/gaidon/progs/dense_trajectory/debug64/FeatTrack"

$exe $video /dev/shm/${video/*\//}.siftgeo $L $W 1 100 $desc $R

ls -lh /dev/shm/$sampid.siftgeo
rm -f /dev/shm/$sampid.siftgeo
