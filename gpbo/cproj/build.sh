#!/usr/bin/env bash

g++ -std=c++11   -c -O3 -fPIC  GPsimple.cpp

g++ -std=c++11   -c -O3 -fPIC  bayesutils.cpp

g++ -std=c++11   -c -O3 -fPIC  direct.cpp

g++ -std=c++11   -c -O3 -fPIC  hypsearch.cpp

g++ -std=c++11   -c -O3 -fPIC  kernel.cpp

g++ -std=c++11   -c -O3 -fPIC  libGP.cpp

g++ -std=c++11   -c -O3 -fPIC  matern.cpp

g++ -std=c++11   -c -O3 -fPIC  misctools.cpp

g++ -std=c++11   -c -O3 -fPIC  newmain.cpp

g++ -std=c++11   -c -O3 -fPIC  simplekernels.cpp

g++ -std=c++11    -o libcproj.so  GPsimple.o bayesutils.o direct.o hypsearch.o kernel.o libGP.o matern.o misctools.o newmain.o simplekernels.o -lblas -llapack -llapacke -shared -fPIC
