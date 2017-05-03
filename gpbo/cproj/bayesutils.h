/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   bayesutils.h
 * Author: mark
 *
 * Created on 15 December 2015, 18:34
 */

#ifndef BAYESUTILS_H
#define BAYESUTILS_H

#include "GPsimple.h"

extern "C" int LCB(GP* g, int n, double* X, int* D, double p, double* R);
extern "C" int EI_gp(GP* g, int n, double* X, int* D, double* R, bool fixI, double II);
extern "C" int lEI_gp(GP* g, int n, double* X, int* D, double* R, bool fixI, double II);

#endif /* BAYESUTILS_H */

