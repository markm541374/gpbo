/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "GPsimple.h"
#include <cmath>
#include <vector>
#include <stdio.h>
#include "misctools.h"

extern "C" int LCB(GP* g, int n, double* X, int* D, double p, double* R){
    std::vector<double> U = std::vector<double>(2*n);
    g->infer_diag(n,X,D,&U[0]);
    for (int i=0; i<n; i++){
        R[i] = U[i] - p*sqrt(U[i+n]);
    }
    return 0;
}

extern "C" int EI_gp(GP* g, int n, double* X, int* D, double* R){
    
    std::vector<double> U = std::vector<double>(2*n);
    g->infer_diag(n,X,D,&U[0]);
    
    for (int i=n; i<2*n; i++){
        U[i] = sqrt(U[i]);
    }
    
    double ymin = 1e99;
    for (int i=0; i<g->N; i++){
        //printf("%f_",g->Yd[i]);
        if(g->Yd[i]<ymin){ymin = g->Yd[i];}
    }
    
    EI(&U[0], &U[n], ymin, n, &R[0]);
    //printf("{%f %f %f %f}",U[0],U[n],R[0],ymin);
    return 0;
}

extern "C" int lEI_gp(GP* g, int n, double* X, int* D, double* R){
    
    std::vector<double> U = std::vector<double>(2*n);
    g->infer_diag(n,X,D,&U[0]);
    
    for (int i=n; i<2*n; i++){
        U[i] = sqrt(U[i]);
    }
    
    double ymin = 1e99;
    for (int i=0; i<g->N; i++){
        //printf("%f_",g->Yd[i]);
        if(g->Yd[i]<ymin){ymin = g->Yd[i];}
    }
    
    lEI(&U[0], &U[n], ymin, n, &R[0]);
    //printf("{%f %f %f %f}",U[0],U[n],R[0],ymin);
    return 0;
}