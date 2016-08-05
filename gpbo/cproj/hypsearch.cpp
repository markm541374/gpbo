/*
 * hypsearch.cpp
 *
 *  Created on: 24 Oct 2015
 *      Author: mark
 */

#include <vector>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include "kernel.h"
#include "direct.h"
#include "GPsimple.h"
#include <cblas.h>
#include <lapacke.h>

const double PI = 3.141592653589793238463;
const double L2PI = log(PI*2);

//GP variables
int D;
int N;
int ki;
int nhyp;
std::vector<double> ih;

std::vector<double> X;
std::vector<double> Y;
std::vector<int> Dx;
std::vector<double> S;

std::vector<double>  Kxx;
std::vector<double> Yd;

std::vector<double> priorm;
std::vector<double> priors;

//default direct hyperparameters

int maxint = 5000;
double fglob = -1e8;

//set direct hyperparameters
extern "C" void SetHypSearchPara(int mx, double fg){
	maxint = mx;
	fglob = fg;
	return;
}
//llk function to call direct from
double llk(int directsearchdim, double* hyp){
        std::vector<double> tmp = std::vector<double>(nhyp);
	int c = hypsearchconvert(hyp,ki,D,&tmp[0]);
        double R;
        bool verbose = false;
        if (verbose){
            printf("At [ ");
            for (int i=0; i<nhyp; i++){printf("%f ",hyp[i]);}
            printf("]: ");
        }
        GP_LKonly G = GP_LKonly(D, N, &X[0], &Y[0], &S[0], &Dx[0], ki, &tmp[0], &R);
        if (verbose){
            printf("%f\n",R);
        }
	return -R;
}

//posterior function to call direct from
double post(int directsearchdim, double* hyp){
    double pst = llk(directsearchdim, hyp);
    
    for (int i=0; i<directsearchdim; i++){
        pst-=0.5*pow((hyp[i]-priorm[i])/priors[i],2);
    }
    return pst;
}
extern "C" int HypSearchMLE(int d, int n, double* Xin, double* Yin, double* Sin, int* Din, double* lb, double* ub, int kernelindex, double* Rhyp, double* lk){

	ki = kernelindex;
        if (ki==2){printf("this linXsquexp has the 4th hyperparameter set to 1 so this is searching an unused dimension");}
	N = n;
	D = d;
        nhyp = numhyp(ki,D);
	ih.resize(D+1);
	X.resize(N*D);
	Y.resize(N);
	S.resize(N);
	Dx.resize(N);
	Kxx.resize(N*N);
	Yd.resize(N);

	for (int i=0; i<N*D; i++){
		X[i] = Xin[i];
	}
	for (int i=0; i<N; i++){
		Y[i] = Yd[i]= Yin[i];
		S[i] = Sin[i];
		Dx[i] = Din[i];
	}
	std::vector<double> xbest = std::vector<double>(nhyp,0.);
        printf("DIRECT searching LLK in %d hyperparameters between (",nhyp);
        for (int i=0; i<nhyp; i++){printf("%f ",lb[i]);}
        printf("\b) and (");
        for (int i=0; i<nhyp; i++){printf("%f ",ub[i]);}
        printf("\b)\n");
        
        direct(nhyp,&lb[0],&ub[0],maxint,fglob,&xbest[0],lk,llk);
        
        
        int c = hypsearchconvert(&xbest[0],ki,D,Rhyp);
        
        printf("DIRECT found (");
        for (int i=0; i<nhyp; i++){printf("%f ",xbest[i]);}
        printf("\b) raw with llk %f which is (",*lk);
        for (int i=0; i<nhyp; i++){printf("%f ",Rhyp[i]);}
        printf("\b)true\n");
        
	return 0;
}

extern "C" int HypSearchMAP(int d, int n, double* Xin, double* Yin, double* Sin, int* Din, double* m, double* s, double z, int kernelindex, double* Rhyp, double* lk){
        //this uses log-normal priors and bounds the search at margin +-z orders of magnitude from the mean
	ki = kernelindex;
        if (ki==2){printf("this linXsquexp has the 4th hyperparameter set to 1 so this is searching an unused dimension");}
	N = n;
	D = d;
        nhyp = numhyp(ki,D);
	ih.resize(D+1);
	X.resize(N*D);
	Y.resize(N);
	S.resize(N);
	Dx.resize(N);
	Kxx.resize(N*N);
	Yd.resize(N);
        priorm.resize(nhyp);
        priors.resize(nhyp);
        
	for (int i=0; i<N*D; i++){
		X[i] = Xin[i];
	}
	for (int i=0; i<N; i++){
		Y[i] = Yd[i]= Yin[i];
		S[i] = Sin[i];
		Dx[i] = Din[i];
	}
        std::vector<double> lbd = std::vector<double>(nhyp);
        std::vector<double> ubd = std::vector<double>(nhyp);
        for (int i=0; i<nhyp; i++){
            priorm[i] = m[i];
            priors[i] = s[i];
            lbd[i] = m[i]-z*s[i];
            ubd[i] = m[i]+z*s[i];                 
        }
	std::vector<double> xbest = std::vector<double>(nhyp,0.);
        printf("DIRECT searching MAP in %d hyperparameters between (",nhyp);
        for (int i=0; i<nhyp; i++){printf("%f ",lbd[i]);}
        printf("\b) and (");
        for (int i=0; i<nhyp; i++){printf("%f ",ubd[i]);}
        printf("\b)\n");
        
        direct(nhyp,&lbd[0],&ubd[0],maxint,fglob,&xbest[0],lk,post);
        
        //printf("zzz");
        int c = hypsearchconvert(&xbest[0],ki,D,Rhyp);
        //printf("zzz");
        printf("DIRECT found (");
        for (int i=0; i<nhyp; i++){printf("%f ",xbest[i]);}
        printf("\b) raw with llk %f which is (",*lk);
        for (int i=0; i<nhyp; i++){printf("%e ",Rhyp[i]);}
        printf("\b)true\n");
        //printf("zzz");
	return 0;
}