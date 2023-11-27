CXXFLAGS = -O3 -Wall -std=gnu++17 -g $(shell pkg-config --cflags fftw3 glib-2.0)

all: designblep testblep

designblep.o: designblep.cc

designblep: designblep.o
	$(CXX) -o designblep designblep.o $(LDFLAGS) $(shell pkg-config --libs fftw3)

testblep.o: testblep.cc

testblep: testblep.o bleposcdata.o
	$(CXX) -o testblep testblep.o bleposcdata.o $(LDFLAGS) $(shell pkg-config --libs fftw3 glib-2.0)

clean:
	rm -f designblep designblep.o testblep.o bleposcdata.o
