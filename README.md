Practical Bayesian Optimization for Variable Cost Objectives
=================================================

This repository contains the implementation of the work in [Practical Bayesian Optimization for Variable Cost Objectives](https://arxiv.org/abs/1703.04335) and [Optimization, fast and slow: optimally switching between local and Bayesian optimization](https://arxiv.org/abs/1805.08610).


Installation
--------------------------------------
The project should install for linux via
```bash
pip install git+https://github.com/markm541374/gpbo@master
```
I reccomend using a virtualenv. On ubuntu 14.04 64bit this should be all that is required. It may be necessary to recompile the code in /gpbo/cproj which is a custom Gaussian Process library allowing observation and inference of first of second derivatives for any kernels that have all the necessary higher order derivatives defined in kernels.ccp (so far only squared exponential and Matern 5/2). 
Assuming g++ is available and all libraries are in the default path this can be done using the script /gpbo/cproj/build.sh . Otherwise the script /gpbo/cproj/buildgrey.sh will need to be modified to suit your system by explicitly providing paths to lapack, cblas and any other libraries that may not be found by defalt.

Examples
--------------------------------------
The /examples folder provides a minimal demonstration of BLOSSOM, PES, PES adapted for variable fidelity objectives, and EI using both slice sampled and maximum likleihood hyperparameters.
