CXXFLAGS=-O3 -Wall -g -std=c++17 -ffast-math

all: abcfit polyfit perffloat testadsr

testadsr.o: testadsr.cc flexadsr.hh

testadsr: testadsr.o
	$(CXX) -o testadsr testadsr.o
