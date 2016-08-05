/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   matern.h
 * Author: mark
 *
 * Created on 29 January 2016, 10:46
 */

#ifndef SIMPLEKERNELS_H
#define SIMPLEKERNELS_H

int mat52perconv(double *h, int D, double* ih);
int mat52persearchconv(double *h, int D, double* ih);
double mat52per(double *x1, double *x2, int d1, int d2, int D, double* ih, double* smodel);

int mat52pptconv(double *h, int D, double* ih);
int mat52pptsearchconv(double *h, int D, double* ih);
double mat52ppt(double *x1, double *x2, int d1, int d2, int D, double* ih, double* smodel);

int matppconv(double *h, int D, double* ih);
int matppsearchconv(double *h, int D, double* ih);
double matpp(double *x1, double *x2, int d1, int d2, int D, double* ih, double* smodel);

int devconv(double *h, int D, double* ih);
int devsearchconv(double *h, int D, double* ih);
double dev(double *x1, double *x2, int d1, int d2, int D, double* ih, double* smodel);
#endif /* SIMPLEKERNELSN_H */

