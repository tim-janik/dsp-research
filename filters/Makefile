CXXFLAGS = -O2 -ffast-math -Wall -g -I /home/stefan/src/pandaresampler/lib

all: testfilter

testfilter.o: testfilter.cc skfilter.hh laddervcf.hh

testfilter: testfilter.o
	$(CXX) -o testfilter testfilter.o $(LDFLAGS)

clean:
	rm -f testfilter testfilter.o
