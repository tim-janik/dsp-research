all: build/DistortionPlugin_artefacts/VST3/DistortionPlugin.vst3/Contents/x86_64-linux/DistortionPlugin.so distortion

build/DistortionPlugin_artefacts/VST3/DistortionPlugin.vst3/Contents/x86_64-linux/DistortionPlugin.so: DistortionPlugin.cpp
	mki.sh

distortion: distortion.cc distortiondsp.hh
	$(CXX) -o distortion distortion.cc -O3 -I/home/stefan/src/pandaresampler/include/ -ftree-vectorize -march=znver1 -mtune=znver1 -ffast-math -g -fno-tree-slp-vectorize -Wall
