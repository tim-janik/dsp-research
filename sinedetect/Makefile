all: sinedetect

sinedetect: sinedetect.cc
	g++ -std=c++17 -Wall -g -O3 -o sinedetect sinedetect.cc `pkg-config --cflags --libs glib-2.0 bse`
