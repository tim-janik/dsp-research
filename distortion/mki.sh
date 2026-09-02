set -e

mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=RelWithAssert
make
cp -avr DistortionPlugin_artefacts/RelWithAssert/VST3/DistortionPlugin.vst3 ~/.vst3
