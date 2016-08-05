/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   newmain.cpp
 * Author: mark
 *
 * Created on 20 June 2016, 11:58
 */

#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cmath>

#include "libGP.h"
#include "kernel.h"
/*
 * 
 */
int main(int argc, char** argv) {
    
    
    int n = 3;
    int Dim = 1;
    printf("data set of dim %d and size %d:\n",Dim,n);
    std::vector<double> x = std::vector<double>({-0.1,0.2,0.5});
    printf("x_data of size %lu:\n[",x.size());
    for (int i=0; i<x.size();i++){
        printf("%f ",x[i]);
    }
    printf("]\n");
    
    std::vector<double> y = std::vector<double>({0.1,0.3,0.2});
    printf("y_data of size %lu:\n[",y.size());
    for (int i=0; i<y.size();i++){
        printf("%f ",y[i]);
    }
    printf("]\n");
    
    std::vector<double> s = std::vector<double>({1e-6,1e-6,1e-6});
    printf("s_data of size %lu:\n[",s.size());
    for (int i=0; i<s.size();i++){
        printf("%f ",s[i]);
    }
    printf("]\n");
    
    std::vector<int> d = std::vector<int>({0,0,0});
    printf("d_data of size %lu:\n[",x.size());
    for (int i=0; i<d.size();i++){
        printf("%d ",d[i]);
    }
    printf("]\n");
    
    std::vector<double> h = std::vector<double>({1.2,0.5});
    std::vector<double> ih = std::vector<double>({1.2,0.5});
    int size = 1;
    printf("h_data of size %lu:\n[",h.size());
    for (int i=0; i<h.size();i++){
        printf("%f ",h[i]);
    }
    printf("]\n");
    
    int ki = 9;
    printf("kernel number: %d\n",ki);
    
    //libGP.hypconvert(self.hyp.ctypes.data_as(ctpd),cint(self.Kindex), cint(self.dim), self.ihyp.ctypes.data_as(ctpd))
    hypconvert(&h[0],ki,Dim,&ih[0]);
    printf("ih_data of size %lu:\n[",ih.size());
    for (int i=0; i<ih.size();i++){
        printf("%f ",ih[i]);
    }
    printf("]\n");
    
    //self.s = libGP.newGP_hypset(cint(self.D),cint(self.n),cint(kf[0].Kindex),X_s.ctypes.data_as(ctpd),Y_s.ctypes.data_as(ctpd),S_s.ctypes.data_as(ctpd),(cint*len(Dx))(*Dx),allhyp.ctypes.data_as(ctpd),cint(self.size))
    int g = newGP_hypset(Dim,n,ki,&x[0],&y[0],&s[0],&d[0],&h[0],size);    
    //ping(g,1);
    presolv(g,1);
    //ping(g,1);
    killGP(g,1);
    //ping(g,1);
    

    return 0;
}

