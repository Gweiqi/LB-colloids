#!/bin/bash

convert -delay 20 -loop 0 *.png synth.gif
ffmpeg -r 30 -f gif -i synth5.gif -codec:v mpeg4 -flags:v +qscale -global_quality:v 0 -codec:a libmp3lame synth5.avi
