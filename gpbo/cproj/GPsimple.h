/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   GPsimple.h
 * Author: mark
 *
 * Created on 26 November 2015, 18:33
 */

#ifndef GPSIMPLE_H
#define GPSIMPLE_H
#include <vector>
class GP_LKonly{
public:
	int D;
	int N;
	int K;
	//double lk;

	GP_LKonly(int d, int n, double* Xin, double* Yin, double* Sin, int* Din, int kindex, double* hyp, double* R);
};


class GP{
public:
    int D;
    int N;
    int K;
    int maxinfer;
    int nhyp;
    std::vector<double> ih;
    std::vector<double>  Kxx;
    std::vector<double> X;
    std::vector<double> Y;
    std::vector<double> Yd;
    std::vector<int> Dx;
    std::vector<double> S;
    std::vector<double> Ksx;
    std::vector<double> Ksx_T;

    GP(int d, int n, int kindex);
    int build_K();
    int set_Y(double* Yin);
    int set_S(double* Sin);
    int set_X(double* Xin);
    int set_D(int* Din);
    int set_hyp(double* hyp);
    void ping();
    int fac();
    int presolv();
    int infer_diag(int Ns, double* Xs, int* Ds, double* R);
    int infer_m(int Ns, double* Xs, int* Ds, double* R);
    int infer_full(int Ns, double* Xs, int* Ds, double* R);
    int draw(int N, double* X, int* D, double*R, int m);
    int llk(double* R);
    int infer_m_partial(int kp, double* ihp, int Ns, double* Xs, int* Ds, double* R);
    int get_cho(double* C);
};
#endif /* GPSIMPLE_H */

