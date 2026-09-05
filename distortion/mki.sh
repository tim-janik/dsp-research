#!/bin/bash
set -e

mkdir -p build
cd build

if [ ! -f CMakeCache.txt ]; then
    cmake .. -DCMAKE_BUILD_TYPE=RelWithAssert
fi

make -j"$(nproc)"
cp -avr DistortionPlugin_artefacts/RelWithAssert/VST3/DistortionPlugin.vst3 ~/.vst3
