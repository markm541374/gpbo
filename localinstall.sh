#!/usr/bin/env bash



(cd gpbo/core; cython *.pyx)
pip install --upgrade --user -e .

mv *.so gpbo/core/