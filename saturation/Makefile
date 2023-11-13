CXXFLAGS = -O3 -ftree-vectorize -march=native -ffast-math -Wall -g -I /home/stefan/src/pandaresampler/lib $(shell pkg-config --cflags fftw3f) -std=c++17
CXX = clang++

all: saturation

saturation.o: saturation.cc saturationdsp.hh

saturation: saturation.o
	$(CXX) -o saturation saturation.o $(LDFLAGS) $(shell pkg-config --libs fftw3f)

clean:
	rm -f saturation saturation.o
