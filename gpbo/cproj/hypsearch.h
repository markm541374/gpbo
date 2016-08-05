/*
 * hypsearch.h
 *
 *  Created on: 24 Oct 2015
 *      Author: mark
 */

#ifndef HYPSEARCH_H_
#define HYPSEARCH_H_

extern "C" int HypSearchMLE(int d, int n, double* Xin, double* Yin, double* Sin, int* Din, double* lb, double* ub, int kernelindex, double* Rhyp, double lk);
extern "C" int HypSearchMAP(int d, int n, double* Xin, double* Yin, double* Sin, int* Din, double* m, double* s, double z, int kernelindex, double* Rhyp, double lk);

extern "C" void SetHypSearchPara(int mx, double fg);


#endif /* HYPSEARCH_H_ */