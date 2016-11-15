#!/usr/bin/env bash

g++ -std=c++11   -c -O2 -fPIC  -MMD -MP -MF "GPsimple.o.d" -o GPsimple.o GPsimple.cpp

g++ -std=c++11   -c -O2 -fPIC  -MMD -MP -MF "bayesutils.o.d" -o bayesutils.o bayesutils.cpp

g++ -std=c++11   -c -O2 -fPIC  -MMD -MP -MF "direct.o.d" -o direct.o direct.cpp

g++ -std=c++11   -c -O2 -fPIC  -MMD -MP -MF "hypsearch.o.d" -o hypsearch.o hypsearch.cpp

g++ -std=c++11   -c -O2 -fPIC  -MMD -MP -MF "kernel.o.d" -o kernel.o kernel.cpp

g++ -std=c++11   -c -O2 -fPIC  -MMD -MP -MF "libGP.o.d" -o libGP.o libGP.cpp

g++ -std=c++11   -c -O2 -fPIC  -MMD -MP -MF "matern.o.d" -o matern.o matern.cpp

g++ -std=c++11   -c -O2 -fPIC  -MMD -MP -MF "misctools.o.d" -o misctools.o misctools.cpp

g++ -std=c++11   -c -O2 -fPIC  -MMD -MP -MF "newmain.o.d" -o newmain.o newmain.cpp

g++ -std=c++11   -c -O2 -fPIC  -MMD -MP -MF "simplekernels.o.d" -o simplekernels.o simplekernels.cpp

g++ -std=c++11    -o libcproj.so GPsimple.o bayesutils.o direct.o hypsearch.o kernel.o libGP.o matern.o misctools.o newmain.o simplekernels.o -lblas -llapack -llapacke -shared -fPIC
