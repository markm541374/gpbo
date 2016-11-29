#include <vector>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include "matern.h"
#include "simplekernels.h"
int cpdec1conv(double *h, int D, double* ih){
    ih[0] = h[0];
    ih[1] = h[1];
    ih[2] = h[2];
    return 0;
}
//origin var and grad var are logspace, offset in unchanged
int cpdec1searchconv(double *h, int D, double* ih){
    ih[0] = pow(10.,h[0]);
    ih[1] = pow(10.,h[1]);
    ih[2] = pow(10.,h[2]);
    return 0;
}
double cpdec1(double *x1, double *x2, int d1, int d2, int D, double* ih, double* smodel){
    //linear kernel in the first axis. ih are v^2 c b^2. Kernel is v^2 *(x1-c)(x2-c)+b^2
    //D is not to be used, only present for consistency in function definition

    if (d1==0 and d2 == 0){
        //no derivatives

        return ih[2]+pow(ih[1],ih[0])/pow(x1[0]+x2[0]+ih[1],ih[0]);
    }
   else{
        return 0.;
    }
}
int dec1conv(double *h, int D, double* ih){
    ih[0] = h[0];
    ih[1] = h[1];
    return 0;
}
//origin var and grad var are logspace, offset in unchanged
int dec1searchconv(double *h, int D, double* ih){
    ih[0] = pow(10.,h[0]);
    ih[1] = pow(10.,h[1]);
    return 0;
}
double dec1(double *x1, double *x2, int d1, int d2, int D, double* ih, double* smodel){
    //linear kernel in the first axis. ih are v^2 c b^2. Kernel is v^2 *(x1-c)(x2-c)+b^2
    //D is not to be used, only present for consistency in function definition

    if (d1==0 and d2 == 0){
        //no derivatives

        return pow(ih[1],ih[0])/pow(x1[0]+x2[0]+ih[1],ih[0]);
    }
   else{
        return 0.;
    }
}

double dec1cs(double *x1, double *x2, int d1, int d2, int D, double* ih, double* smodel){

    smodel[0] = ih[2];
    return dec1(x1,x2,d1,d2,D,ih,smodel);
}

int dec1csconv(double *h, int D, double* ih){
    ih[0] = pow(10.,h[0]);
    ih[1] = pow(10.,h[1]);
    ih[2] = pow(10.,h[2]);
    return 0;
}
//all are searched in log space
int dec1cssearchconv(double *h, int D, double* ih){
    ih[0] = pow(10.,h[0]);
    ih[1] = pow(10.,h[1]);
    ih[2] = pow(10.,h[2]);
    return 0;
}

int lin1conv(double *h, int D, double* ih){
    ih[0] = pow(h[0],2);
    ih[1] = h[1];
    ih[2] = pow(h[2],2);
    return 0;
}
//origin var and grad var are logspace, offset in unchanged
int lin1searchconv(double *h, int D, double* ih){
    ih[0] = pow(10.,h[0]);
    ih[1] = h[1];
    ih[2] = pow(10.,h[2]);
    return 0;
}
double lin1(double *x1, double *x2, int d1, int d2, int D, double* ih, double* smodel){
    //linear kernel in the first axis. ih are v^2 c b^2. Kernel is v^2 *(x1-c)(x2-c)+b^2
    //D is not to be used, only present for consistency in function definition
    
    if (d1==0 and d2 == 0){
        //no derivatives
        
        return ih[0]*(x1[0]-ih[1])*(x2[0]-ih[1])+ih[2];
    }
    else if (d1==1 and d2==0){
        return ih[0]*(x2[0]-ih[1]);
    }
    else if (d1==0 and d2==1){
        return ih[0]*(x1[0]-ih[1]);
    }
    else if (d1==1 and d2 ==1){
        return ih[0];
    }
    else{
        return 0.;
    }
}

//output scale is squared, lengths are 1/square
int squexpconv(double *h, int D, double* ih){
    ih[0] = pow(h[0],2);
    for (int i=1; i<D+1; i++){
	ih[i] = 1./pow(h[i],2);
    }
    return 0;
}

//all are searched in log space
int squexpsearchconv(double *h, int D, double* ih){
    
    for (int i=0; i<D+1; i++){
        //printf("_%f",h[i]);
	ih[i] = pow(10.,h[i]);
    }
    return 0;
}

double squexp(double *x1, double *x2, int d1, int d2, int D, double* ih, double* smodel){
	double expon = 0.;
	for (int i=0; i<D; i++){
		expon-=pow((x1[i]-x2[i]),2)*ih[i+1];
	}
	double core = ih[0]*exp(0.5*expon);
	if (d1==0 and d2==0){
		//no derivatives
		return core;
	}
        //else{printf("%d %d",d1,d2);}
	std::vector<int> V = std::vector<int>(D,0);
	div_t v1;
	div_t v2;
	int S = 0;

	int sign = 1;
	for (int i=0; i<D; i++){
		v1 = div(d1,pow(8,D-i-1));
		V[D-i-1] += v1.quot;
		S+=v1.quot;
		d1 = v1.rem;
		v2 = div(d2,pow(8,D-i-1));
		V[D-i-1] += v2.quot;
		sign += v2.quot;
		S+=v2.quot;
		d2 = v2.rem;
	}
	int P = 1;
	for (int i=0; i<D; i++){
		P*=V[i]+1;
	}
	sign = 2*(sign%2) -1;
	if (S==1){
		//first derivative
		int i =0;
		while (V[i]==0){i+=1;}
		return -ih[i+1]*(x1[i]-x2[i])*double(sign)*core;
	}
	else if (S==2){
		if (P==3){
			//second derivative
			int i =0;
			while (V[i]==0){i+=1;}
			return ih[i+1]*(ih[i+1]*pow((x1[i]-x2[i]),2)-1.) * double(sign)*core;
		}
		else if (P==4){
			//first derivative on two axes
			int i = 0;
			while (V[i]==0){i+=1;}
			int j = i+1;
			while (V[j]==0){j+=1;}
			//printf("%d,%d\n",i,j);
			return ih[j+1]*(x1[j]-x2[j])*ih[i+1]*(x1[i]-x2[i])*double(sign)*core;
		}
		else{
			printf("invalid derivatives %d %d",d1,d2);
			return 0.;
		}
	}
	else if (S==3){
		if (P==4){
			//third derivative
			//(li**2)*(3*xi-li*xi**3)
			int i = 0;
			while (V[i]==0){i+=1;}
			double x = (x1[i]-x2[i]);
			double l = ih[i+1];
			return pow(l,2)*(3.*x-l*pow(x,3))*double(sign)*core;
		}
		else if (P==6){
			//second and first derivative
			//-xi*li*lj*(lj*xj**2-1)
			//j is the repeated axis
			int i = 0;
			while (V[i]!=1){i+=1;}
			int j = 0;
			while (V[j]!=2){j+=1;}
			//printf("[%d,%d]\n",i,j);
			double xi = (x1[i]-x2[i]);
			double li = ih[i+1];
			double xj = (x1[j]-x2[j]);
			double lj = ih[j+1];
			return -xi*li*lj*(lj*pow(xj,2)-1)*double(sign)*core;
		}
		else if (P==8){
			//three first derivatives
			int i = 0;
			while (V[i]==0){i+=1;}
			int j = i+1;
			while (V[j]==0){j+=1;}
			int k = j+1;
			while (V[k]==0){k+=1;}
			//printf("%d,%d,%d\n",i,j,k);
			//- xi*li*xj*lj*xk*lk
			return -ih[k+1]*(x1[k]-x2[k])*ih[j+1]*(x1[j]-x2[j])*ih[i+1]*(x1[i]-x2[i])*double(sign)*core;

			return 0.;
		}
		else{
			printf("invalid derivatives %d %d",d1,d2);
			return 0.;
		}
	}
	else if(S==4){
		if (P==16){
			//four first derivatives
			//xi*li*xj*lj*xk*lk*xl*ll
			int i = 0;
			while (V[i]==0){i+=1;}
			int j = i+1;
			while (V[j]==0){j+=1;}
			int k = j+1;
			while (V[k]==0){k+=1;}
			int l = k+1;
			while (V[k]==0){l+=1;}
			return ih[l+1]*(x1[l]-x2[l])*ih[k+1]*(x1[k]-x2[k])*ih[j+1]*(x1[j]-x2[j])*ih[i+1]*(x1[i]-x2[i])*double(sign)*core;
		}
		else if (P==12){
			//one second and two first
			//lk*(lk*xk**2-1)*xi*li*xj*lj
			//k is the repeated axis

			int i = 0;
			while (V[i]!=1){i+=1;}
			int j = i+1;
			while (V[j]==0){j+=1;}
			int k = 0;
			while (V[k]!=2){k+=1;}
			double xi = (x1[i]-x2[i]);
			double li = ih[i+1];
			double xj = (x1[j]-x2[j]);
			double lj = ih[j+1];
			double xk = (x1[k]-x2[k]);
			double lk = ih[k+1];

			return lk*(lk*pow(xk,2)-1)*xi*li*xj*lj*double(sign)*core;
		}
		else if (P==9){
			//two second derivatives
			//li*(li*xi**2-1)*lj*(lj*xj**2-1)
			int i = 0;
			while (V[i]==0){i+=1;}
			int j = i+1;
			while (V[j]==0){j+=1;}
			double xi = (x1[i]-x2[i]);
			double li = ih[i+1];
			double xj = (x1[j]-x2[j]);
			double lj = ih[j+1];

			return li*(li*pow(xi,2)-1)*lj*(lj*pow(xj,2)-1)*double(sign)*core;
		}
		else if (P==8){
			//third and first derivative
			//-li*xi*(lj**2)*(3*xj-lj*xj**3)
			//j is the repeated axis
			int i = 0;
			while (V[i]!=1){i+=1;}
			int j = 0;
			while (V[j]!=3){j+=1;}
			double xi = (x1[i]-x2[i]);
			double li = ih[i+1];
			double xj = (x1[j]-x2[j]);
			double lj = ih[j+1];
			return -li*xi*pow(lj,2)*(3*xj-lj*pow(xj,3))*double(sign)*core;
		}
		else if (P==5){
			//fourth derivative
			//(li*xi)**4 - 6*(li**3)*(xi**2) + 3*li**2
			int i = 0;
			while (V[i]==0){i+=1;}
			double xi = (x1[i]-x2[i]);
			double li = ih[i+1];
			return (pow(li*xi,4) - 6*pow(li,3)*pow(xi,2) + 3*pow(li,2))*double(sign)*core;
		}
		else{
			printf("invalid derivatives %d %d",d1,d2);
			return 0.;
		}
	}
	else{
		printf("invalid derivativesx %d %d %d",d1,d2,S);
		return 0.;
	}

}
/*
 --------------------------------------------------------------------------
 
 */
double squexpcs(double *x1, double *x2, int d1, int d2, int D, double* ih, double* smodel){
    
    smodel[0] = ih[D+1];
    return squexp(x1,x2,d1,d2,D,ih,smodel);
}
            
int squexpcsconv(double *h, int D, double* ih){
    ih[0] = pow(h[0],2);
    for (int i=1; i<D+1; i++){
	ih[i] = 1./pow(h[i],2);
        
    }
    ih[D+1] = h[D+1];
    return 0;
}
//all are searched in log space
int squexpcssearchconv(double *h, int D, double* ih){
    for (int i=0; i<D+2; i++){
        //printf("_%f",h[i]);
	ih[i] = pow(10.,h[i]);
    }
    
    return 0;
}
/*
 --------------------------------------------------------------------------
 
 */
double squexpps(double *x1, double *x2, int d1, int d2, int D, double* ih, double* smodel){
    //noise is parabolic in the first axis
    //I use x1 as s is only used when s1=s2
    double tmp = -ih[D+1]*(x1[0]-ih[D+2])*(x1[0]-ih[D+3]);
    if (tmp>0.){smodel[0] = tmp;}
    else {smodel[0]=1e-9;}
    
    return squexp(x1,x2,d1,d2,D,ih,smodel);
}
int squexppsconv(double *h, int D, double* ih){
    ih[0] = pow(h[0],2);
    for (int i=1; i<D+1; i++){
	ih[i] = 1./pow(h[i],2);
        
    }
    ih[D+1] = h[D+1];
    ih[D+2] = h[D+2];
    ih[D+3] = h[D+3];    
    return 0;
}
//all are searched in log space exept last two which are x coordinates
int squexppssearchconv(double *h, int D, double* ih){
    for (int i=0; i<D+2; i++){
        //printf("_%f",h[i]);
	ih[i] = pow(10.,h[i]);
    }
    ih[D+2] = h[D+2];
    ih[D+3] = h[D+3];
    return 0;
}
/*
 --------------------------------------------------------------------------
 
 */
double squexpbs(double *x1, double *x2, int d1, int d2, int D, double* ih, double* smodel){
    //noise is parabolic in the first axis
    //I use x1 as s is only used when s1=s2
    double tmp = ih[D+1]*pow(x1[0],ih[D+2]*ih[D+3])*pow(1.-x1[0],ih[D+2]*(1.-ih[D+3]))+1e-10;
    if (tmp>0.){smodel[0] = tmp;}
    else {smodel[0]=1e-9;}
    
    return squexp(x1,x2,d1,d2,D,ih,smodel);
}
int squexpbsconv(double *h, int D, double* ih){
    ih[0] = pow(h[0],2);
    for (int i=1; i<D+1; i++){
	ih[i] = 1./pow(h[i],2);
        
    }
    ih[D+1] = h[D+1];
    ih[D+2] = h[D+2];
    ih[D+3] = h[D+3];    
    return 0;
}

//all are searched in log space exept last one which is ~sigmoid
int squexpbssearchconv(double *h, int D, double* ih){
    for (int i=0; i<D+3; i++){
        //printf("_%f",h[i]);
	ih[i] = pow(10.,h[i]);
    }
    ih[D+3] = 1./(1.+pow(2.718,h[D+3]));
    return 0;
}
/*
 --------------------------------------------------------------------------
 
 */
int linXPsquexpsearchconv(double *h, int D, double* ih){
    lin1searchconv(h,D,ih);
    squexpsearchconv(&h[3],D-1,&ih[3]);
    ih[3]=1.;
    return 0;
}

int linXPsquexpconv(double *h, int D, double* ih){
    lin1conv(h,D,ih);
    squexpconv(&h[3],D-1,&ih[3]);
    ih[3]=1.;
    return 0;
}
double linXPsquexp(double *x1, double *x2, int d1, int d2, int D, double* ih, double* smodel){
    //linear kernel multiplied by squared exp in perpendicular directions
    //linear in the first axis, but with varying m and c in the other dimensions
    //linear hyps are first 3 hyperparameters
    //I can only not do product rule because these are entirely in different axes
    return lin1(x1,x2,d1,d2,D,ih,smodel) * squexp(&x1[1],&x2[1],d1/8,d2/8,D-1,&ih[3],smodel);
    
}
int linsquexpXPsquexpsearchconv(double *h, int D, double* ih){
    lin1searchconv(h,D,ih);
    squexpsearchconv(&h[3],1,&ih[3]);
    ih[3]=1.;
    squexpsearchconv(&h[3],D-1,&ih[5]);
    ih[5]=1.;
    return 0;
}
int linsquexpXPsquexpconv(double *h, int D, double* ih){
    lin1conv(h,D,ih);
    squexpconv(&h[3],1,&ih[3]);
    ih[3]=1.;
    squexpconv(&h[5],D-1,&ih[5]);
    ih[5]=1.;
    return 0;
}

double linsquexpXPsquexp(double *x1, double *x2, int d1, int d2, int D, double* ih, double* smodel){
    //if (d1!=0 or d2!=0){printf("xxx%d %dxxx",d1,d2);}
    return (lin1(x1,x2,d1,d2,D,ih,smodel) + squexp(x1,x2,d1%8,d2%8,1,&ih[3],smodel))* squexp(&x1[1],&x2[1],d1/8,d2/8,D-1,&ih[5],smodel);
    
 
}

int squexp1Ssquexpconv(double *h, int D, double* ih){
    squexpconv(h, D, ih);
    squexpconv(&h[D+1], 1, &ih[D+1]);
    return 0;
}

//all are searched in log space
int squexp1Ssquexpsearchconv(double *h, int D, double* ih){
    squexpsearchconv(h, D, ih);
    squexpsearchconv(&h[D+1], 1, &ih[D+1]);
    return 0;
}
//squexp in D dimensions plus an additional
double squexp1Ssquexp(double *x1, double *x2, int d1, int d2, int D, double* ih, double* smodel){
    return squexp(x1, x2, d1, d2, D, ih,smodel)+ squexp(x1,  x2, d1%8, d2%8, 1, &ih[D+1],smodel);
}

int squexpsquexpPsquexpconv(double *h, int D, double* ih){
    squexpconv(h, 1, ih);
    squexpconv(&h[2], D, &ih[2]);
    
    return 0;
}

//all are searched in log space
int squexpsquexpPsquexpsearchconv(double *h, int D, double* ih){
    squexpsearchconv(h, 1, ih);
    squexpsearchconv(&h[2], D, &ih[2]);
    
    
    return 0;
}
//squexp +squexp in 1st dimension * squexp in teh others
double squexpsquexpPsquexp(double *x1, double *x2, int d1, int d2, int D, double* ih, double* smodel){
    std::vector<double> tmp = std::vector<double>(D);
    tmp[0]=1.;
    for (int i=0; i<D-1;i++){tmp[i+1]=ih[4+i];}
    
    return (squexp(x1, x2, d1%8, d2%8, 1, ih, smodel)+squexp(x1, x2, d1%8, d2%8, 1, &ih[2],smodel))*squexp(&x1[1], &x2[1], d1/8, d2/8, D-1, &tmp[0],smodel) ;
}


typedef double (*KP)(double*, double*, int, int, int, double*,double*);

extern "C" const KP kern[18] = {&squexp,&lin1,&linXPsquexp,&linsquexpXPsquexp,squexp1Ssquexp,&squexpsquexpPsquexp,&squexpcs,&squexpps,&squexpbs,&mat52,&mat52cs,&mat52per,&mat52ppt,&dev,&matpp,&dec1,&dec1cs,&cpdec1};

extern "C" double k(double *x1, double *x2, int d1, int d2, int D, double* ih, int kindex, double* smodel){
    smodel[0] = 0.;
return kern[kindex](&x1[0], &x2[0], d1, d2, D, &ih[0],smodel);
}
//double (*kern)(double *x1, double *x2, int d1, int d2, int D, double* ih) = &k;



typedef int (*HP)(double *h, int D, double* ih);

extern "C" const HP hypcons[18] = {&squexpconv,&lin1conv,&linXPsquexpconv,&linsquexpXPsquexpconv,&squexp1Ssquexpconv,&squexpsquexpPsquexpconv,&squexpcsconv,&squexppsconv,&squexpbsconv,&mat52conv,&mat52csconv,&mat52perconv,&mat52pptconv,&devconv,&matppconv,&dec1conv,&dec1csconv,&cpdec1conv};

extern "C" int hypconvert(double *h, int kindex, int D, double* ih){
    return hypcons[kindex](h,D,ih);
}

typedef int (*SP)(double *h, int D, double* ih);

extern "C" const SP hypsearchcons[18] = {&squexpsearchconv,&lin1searchconv,&linXPsquexpsearchconv,&linsquexpXPsquexpsearchconv,&squexp1Ssquexpsearchconv,&squexpsquexpPsquexpsearchconv,&squexpcssearchconv,&squexppssearchconv,&squexpbssearchconv,&mat52searchconv,&mat52cssearchconv,&mat52persearchconv,&mat52pptsearchconv,&devsearchconv,&matppsearchconv,&dec1searchconv,&dec1cssearchconv,&cpdec1searchconv};

extern "C" int hypsearchconvert(double *h, int kindex, int D, double* ih){
    return hypsearchcons[kindex](h,D,ih);
}
extern "C" int numhyp(int kindex, int D){
    if (kindex==0){
        return D+1;
    }
    else if (kindex==1){
        return 3;
    }
    else if (kindex==2){
        return D+3;
    }
    else if (kindex==3){
        return D+5;
    }
    else if (kindex==4){
        return D+3;
    }
    else if (kindex==5){
        return D+3;
    }
    else if (kindex==6){
        return D+2;
    }
    else if (kindex==7){
        return D+4;
    }
    else if (kindex==8){
        return D+4;
    }
    else if (kindex==9){
        return D+1;
    }
    else if (kindex==10){
        return D+2;
    }
    else if (kindex==11){
        return 3;
    }
    else if (kindex==12){
        return 4;
    }
    else if (kindex==13){
        //dev
        return 8;
    }
    else if (kindex==14){
        //dev
        return 5;
    }
    else if (kindex==15){
        return 2;
    }
    else if (kindex==16){
        return 3;
    }
    else if (kindex==17){
        return 3;
    }
    else{
        printf("%d %d a bad thing happened :(",kindex,D);
        return -1;
    }
}