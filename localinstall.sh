#!/usr/bin/env bash

cp gpbo/cproj/dist/Release/GNU-Linux/libcproj.so gpbo/cproj/

(cd gpbo/core; cython *.pyx)
pip install --upgrade --user -e .

mv *.so gpbo/core/