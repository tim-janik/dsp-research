CXXFLAGS = -O2 -ffast-math -Wall -g -I /home/stefan/src/pandaresampler/lib $(shell pkg-config --cflags fftw3f)

all: testfilter testtanh

testfilter.o: testfilter.cc skfilter.hh laddervcf.hh

testtanh.o: testtanh.cc

testfilter: testfilter.o
	$(CXX) -o testfilter testfilter.o $(LDFLAGS) $(shell pkg-config --libs fftw3f)

testtanh: testtanh.o
	$(CXX) -o testtanh testtanh.o $(LDFLAGS) $(shell pkg-config --libs fftw3f)
clean:
	rm -f testfilter testfilter.o testtanh testtanh.o
