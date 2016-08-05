/*
 * kernel.h
 *
 *  Created on: 24 Oct 2015
 *      Author: mark
 */

#ifndef KERNEL_H_
#define KERNEL_H_
const int n=15; //how many kernels there are


//k is the kindex^th kernel in d dimensions between x1 and x2 and differentiated accoding to d1 and d2 with adjusted hyperparameters ih
//d1/d2 denote derivatives of the kernel. the remainder of d/8^i is the number of derivativees in the i^th axis
extern "C" double k(double *x1, double *x2, int d1, int d2, int D, double* ih, int kindex, double* smodel);
//extern "C" double (*kern)(double *x1, double *x2, int d1, int d2, int D, double* ih);
typedef double (*KP)(double *x1, double *x2, int d1, int d2, int D, double* ih, double* smodel);
extern "C" const KP kern[n];


//convert the base hyperparameters to form used by kernel function (eg: lengthscale becomes 1/lenghtscale^2)
extern "C" int hypconvert(double *h, int kindex, int D, double* ih);
typedef int (*HP)(double *h, int D, double* ih);
extern "C" const HP hypcons[n];

//convert transform from hyps used in search to base form (eg: lengthscales are searched in log space)
extern "C" int hypsearchconvert(double *h, int kindex, int D, double* ih);
typedef int (*SP)(double *h, int D, double* ih);
extern "C" const SP hypcons[n];


extern "C" int numhyp(int kindex, int D);
#endif /* KERNEL_H_ */