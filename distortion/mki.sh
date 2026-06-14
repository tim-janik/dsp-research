set -e

mkdir -p build
cd build
cmake ..
make
cp -avr DistortionPlugin_artefacts/VST3/DistortionPlugin.vst3 ~/.vst3
