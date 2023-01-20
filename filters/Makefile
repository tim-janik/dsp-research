CXXFLAGS = -O2 -ffast-math -Wall -g -I /home/stefan/src/pandaresampler/lib $(shell pkg-config --cflags fftw3f)

all: testfilter

testfilter.o: testfilter.cc skfilter.hh laddervcf.hh

testfilter: testfilter.o
	$(CXX) -o testfilter testfilter.o $(LDFLAGS) $(shell pkg-config --libs fftw3f)

clean:
	rm -f testfilter testfilter.o
