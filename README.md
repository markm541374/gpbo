Practical Bayesian Optimization for Variable Cost Objectives
==================================================

This repository contains the implementation of the work in [Practical Bayesian Optimization for Variable Cost Objectives](https://arxiv.org/abs/1703.04335)


Installation
--------------------------------------
The project should install for linux via
```bash
pip install git+https://github.com/markm541374/gpbo@master
```
I reccomend using a virtualenv. On ubuntu 14.04 64bit this should be all that is required. It may be necessary to recompile the code in /gpbo/cproj . 
Assuming g++ is available and all libraries are in the default path this can be done using the script /gpbo/cproj/build.sh . Otherwise the script /gpbo/cproj/buildgrey.sh will need to be modified to suit your system by explicitly providing paths to lapack, cblas and any other libraries that may not be found by defalt.

Examples
--------------------------------------
The /examples folder provides a minimal demonstration of the variable fidelity optimization plus PES and expected improvement with both maximum likelihood and samples hyperparameters
