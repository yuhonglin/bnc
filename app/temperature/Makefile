all: sample diagnostic

sample: temperature.cpp
	g++ -O3 -fopenmp -I../../src/ -I/home/honglin/local/lib/R/include/ -L/home/honglin/local/lib/R/lib/ -L/home/honglin/local/lib/ -std=c++11 temperature.cpp /home/honglin/Source/R-3.2.5/src/main/libR.a ../../src/dist/brmath/build/libbrmath.a  ../../src/rng/BasicRNG.cc ../../src/util/logger.cc -o $@ -lRblas -lRlapack -ldl -lreadline -llzma -lpcre -lpng -ljpeg -lz -lpthread


diagnostic: diagnostic.cpp
	g++ -O3 -I../../src/ -I/opt/local/Library/Frameworks/R.framework/Headers/ -L/opt/local/Library/Frameworks/R.framework/Libraries/ -std=c++11 diagnostic.cpp ../../src/dist/brmath/build/libbrmath.a ../../src/rng/BasicRNG.cc ../../src/util/logger.cc -lR -o diagnostic
