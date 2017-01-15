#!/usr/bin/env bash

g++ -std=c++11   -c -O2 -fPIC  -MMD -MP -MF "GPsimple.o.d" -o GPsimple.o GPsimple.cpp -I/usr/include -I/usr/include/openblas -I/usr/include/lapacke


g++ -std=c++11   -c -O2 -fPIC  -MMD -MP -MF "bayesutils.o.d" -o bayesutils.o bayesutils.cpp -I/usr/include -I/usr/include/openblas -I/usr/include/lapacke

g++ -std=c++11   -c -O2 -fPIC  -MMD -MP -MF "direct.o.d" -o direct.o direct.cpp -I/usr/include -I/usr/include/openblas -I/usr/include/lapacke

g++ -std=c++11   -c -O2 -fPIC  -MMD -MP -MF "hypsearch.o.d" -o hypsearch.o hypsearch.cpp -I/usr/include -I/usr/include/openblas -I/usr/include/lapacke


g++ -std=c++11   -c -O2 -fPIC  -MMD -MP -MF "kernel.o.d" -o kernel.o kernel.cpp -I/usr/include -I/usr/include/openblas -I/usr/include/lapacke

g++ -std=c++11   -c -O2 -fPIC  -MMD -MP -MF "libGP.o.d" -o libGP.o libGP.cpp -I/usr/include -I/usr/include/openblas -I/usr/include/lapacke


g++ -std=c++11   -c -O2 -fPIC  -MMD -MP -MF "matern.o.d" -o matern.o matern.cpp -I/usr/include -I/usr/include/openblas -I/usr/include/lapacke

g++ -std=c++11   -c -O2 -fPIC  -MMD -MP -MF "misctools.o.d" -o misctools.o misctools.cpp -I/usr/include -I/usr/include/openblas -I/usr/include/lapacke


g++ -std=c++11   -c -O2 -fPIC  -MMD -MP -MF "newmain.o.d" -o newmain.o newmain.cpp -I/usr/include -I/usr/include/openblas -I/usr/include/lapacke

g++ -std=c++11   -c -O2 -fPIC  -MMD -MP -MF "simplekernels.o.d" -o simplekernels.o simplekernels.cpp -I/usr/include -I/usr/include/openblas -I/usr/include/lapacke

g++ -std=c++11    -o libcproj.so GPsimple.o bayesutils.o direct.o hypsearch.o kernel.o libGP.o matern.o misctools.o newmain.o simplekernels.o -L/usr/lib64  -l:liblapacke.so -l:liblapack.so -l:libgslcblas.so -l:libopenblas.so -lgfortran -shared -fPIC

