/*
 * 1d undifferentiated kernels
 * 
 * 
 */
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>

const double PI = 3.14159265359;
const double OVERPI = 0.31830988618;
const double FIVETHIRDS = 5./3.;
const double SQRT5 = sqrt(5.);
int mat52perconv(double *h, int D, double* ih){
    ih[0] = pow(h[0],2);
    ih[1] = 1./h[1];
    ih[2] = h[2];
    return 0;
}

//all are searched in log space
int mat52persearchconv(double *h, int D, double* ih){
    
    ih[0] = pow(10.,h[0]);
    ih[1] = pow(10.,h[1]);
    ih[2] = pow(10.,h[2]); 
    return 0;
}

double mat52per(double *x1, double *x2, int d1, int d2, int D, double* ih, double* smodel){
        //strictly 1d kernel
        //A mat52( sin2(2pix/p) )
	//ihyps outputscale**2,1/lengthscale,period
	if (d1!=0 or d2!=0 or D!=1){
		printf("no derivatives for this kernel, D must be 1 %d %d %d",d1,d2,D);
		return 0.;
	}
        
	//double r = abs(ih[2]*ih[1]*OVERPI*sin(PI*(x1[0]-x2[0])/ih[2]));
        double r = fabs(double(ih[2]*ih[1]*OVERPI*sin(PI*(x1[0]-x2[0])/ih[2])));
        
        return ih[0]*(1.+SQRT5*r+FIVETHIRDS*pow(r,2))*exp(-SQRT5*r);

}

int mat52pptconv(double *h, int D, double* ih){
    ih[0] = pow(h[0],2);
    ih[1] = 1./h[1];
    ih[2] = h[2];
    ih[3] = 1./pow(h[3],2);
    return 0;
}

//all are searched in log space
int mat52pptsearchconv(double *h, int D, double* ih){
    
    ih[0] = pow(10.,h[0]);//mat52( sin2(2pix/p) )
    ih[1] = pow(10.,h[1]);
    ih[2] = pow(10.,h[2]);
    ih[3] = pow(10.,h[3]); 
    return 0;
}

double mat52ppt(double *x1, double *x2, int d1, int d2, int D, double* ih, double* smodel){
        //strictly 1d kernel
        //A mat52( sin2(2pix/p) )*squexp
	//ihyps outputscale**2,1/lengthscale,period,1/decay**2
	if (d1!=0 or d2!=0 or D!=1){
		printf("no derivatives for this kernel, D must be 1 %d %d %d",d1,d2,D);
		return 0.;
	}
        double dx = x1[0]-x2[0];
	//double r = abs(ih[2]*ih[1]*OVERPI*sin(PI*(x1[0]-x2[0])/ih[2]));
        double r = fabs(double(ih[2]*ih[1]*OVERPI*sin(PI*dx/ih[2])));
        
        return ih[0]*(1.+SQRT5*r+FIVETHIRDS*pow(r,2))*exp(-SQRT5*r)*exp(-0.5*pow(dx,2)*ih[3]);

}

int matppconv(double *h, int D, double* ih){
    ih[0] = pow(h[0],2);
    ih[1] = 1./h[1];
    ih[2] = h[2];
    
    ih[3] = pow(h[3],2);
    ih[4] = 1./h[4];
    
    //ih[6] = pow(h[6],2);
    //ih[7] = 1./h[7];
    
    return 0;
}

//all are searched in log space
int matppsearchconv(double *h, int D, double* ih){
    
    ih[0] = pow(10.,h[0]);
    ih[1] = pow(10.,h[1]);
    ih[2] = pow(10.,h[2]);
    ih[3] = pow(10.,h[3]);
    ih[4] = pow(10.,h[4]);
    
    return 0;
}

double matpp(double *x1, double *x2, int d1, int d2, int D, double* ih, double* smodel){
        //strictly 1d kernel
        //A mat52( sin2(2pix/p) )*squexp
	//ihyps outputscale**2,1/lengthscale,period,1/decay**2
	if (d1!=0 or d2!=0 or D!=1){
		printf("no derivatives for this kernel, D must be 1 %d %d %d",d1,d2,D);
		return 0.;
	}
        double dx = x1[0]-x2[0];
	
        double r1 = fabs(double(ih[2]*ih[1]*OVERPI*sin(PI*dx/ih[2])));
        double mp1 = ih[0]*(1.+SQRT5*r1+FIVETHIRDS*pow(r1,2))*exp(-SQRT5*r1);
        
        
        double r3 = fabs(dx*ih[4]);
        double mm = ih[3]*(1.+SQRT5*r3+FIVETHIRDS*pow(r3,2))*exp(-SQRT5*r3);
        
        
        return mp1+mm;

}

int devconv(double *h, int D, double* ih){
    ih[0] = pow(h[0],2);
    ih[1] = 1./h[1];
    ih[2] = h[2];
    ih[3] = pow(h[3],2);
    ih[4] = 1./h[4];
    ih[5] = h[5];
    ih[6] = pow(h[6],2);
    ih[7] = 1./h[7];
    
    //ih[6] = pow(h[6],2);
    //ih[7] = 1./h[7];
    
    return 0;
}

//all are searched in log space
int devsearchconv(double *h, int D, double* ih){
    
    ih[0] = pow(10.,h[0]);
    ih[1] = pow(10.,h[1]);
    ih[2] = pow(10.,h[2]);
    ih[3] = pow(10.,h[3]);
    ih[4] = pow(10.,h[4]);
    ih[5] = pow(10.,h[5]);
    ih[6] = pow(10.,h[4]);
    ih[7] = pow(10.,h[5]);
    return 0;
}

double dev(double *x1, double *x2, int d1, int d2, int D, double* ih, double* smodel){
        //strictly 1d kernel
        //A mat52( sin2(2pix/p) )*squexp
	//ihyps outputscale**2,1/lengthscale,period,1/decay**2
	if (d1!=0 or d2!=0 or D!=1){
		printf("no derivatives for this kernel, D must be 1 %d %d %d",d1,d2,D);
		return 0.;
	}
        double dx = x1[0]-x2[0];
	
        double r1 = fabs(double(ih[2]*ih[1]*OVERPI*sin(PI*dx/ih[2])));
        double mp1 = ih[0]*(1.+SQRT5*r1+FIVETHIRDS*pow(r1,2))*exp(-SQRT5*r1);
        
        double r2 = fabs(double(ih[5]*ih[4]*OVERPI*sin(PI*dx/ih[5])));
        double mp2 = ih[3]*(1.+SQRT5*r2+FIVETHIRDS*pow(r2,2))*exp(-SQRT5*r2);
        
        double r3 = fabs(dx*ih[7]);
        double mm = ih[6]*(1.+SQRT5*r3+FIVETHIRDS*pow(r3,2))*exp(-SQRT5*r3);
        
        
        return mp1+mp2+mm;

}