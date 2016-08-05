/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   misctools.h
 * Author: markm
 *
 * Created on 26 November 2015, 14:50
 */

#ifndef MISCTOOLS_H
#define MISCTOOLS_H

extern "C" int EI(double*, double*, double, int, double*);
extern "C" int lEI(double*, double*, double, int, double*);
double pdf(double);
double lpdf(double);
double cdf(double);
extern "C" int drawcov(double* K, int n, double* R, int m);
extern "C" int drawk(double* K, int n, double* R, int m);
#endif /* MISCTOOLS_H */

