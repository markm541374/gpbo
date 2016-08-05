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

#ifndef MATERN_H
#define MATERN_H

int mat52conv(double *h, int D, double* ih);
int mat52searchconv(double *h, int D, double* ih);
double mat52(double *x1, double *x2, int d1, int d2, int D, double* ih, double* smodel);
int mat52csconv(double *h, int D, double* ih);
int mat52cssearchconv(double *h, int D, double* ih);
double mat52cs(double *x1, double *x2, int d1, int d2, int D, double* ih, double* smodel);

#endif /* MATERN_H */

