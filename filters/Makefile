CXXFLAGS = -O2 -ffast-math -Wall -g -I /home/stefan/src/pandaresampler/lib $(shell pkg-config --cflags fftw3f)
LDLIBS = $(shell pkg-config --libs fftw3f)

all: testfilter

testfilter: testfilter.cc
