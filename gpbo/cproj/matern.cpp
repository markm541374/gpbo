/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
const double FIVETHIRDS = 5./3.;
const double SQRT5 = sqrt(5.);
int mat52conv(double *h, int D, double* ih){
    ih[0] = pow(h[0],2);
    for (int i=1; i<D+1; i++){
	ih[i] = 1./pow(h[i],2);
    }
    return 0;
}

//all are searched in log space
int mat52searchconv(double *h, int D, double* ih){
    
    for (int i=0; i<D+1; i++){
        //printf("_%f",h[i]);
	ih[i] = pow(10.,h[i]);
    }
    return 0;
}

double mat52(double *x1, double *x2, int d1, int d2, int D, double* ih, double* smodel){
	double r2 = 0.;
	for (int i=0; i<D; i++){
		r2+=pow((x1[i]-x2[i]),2)*ih[i+1];
	}
        double sq5r2 = SQRT5*sqrt(r2);
        double oversq5r2 = 1./sq5r2;
        double core = exp(-sq5r2);
	
	if (d1==0 and d2==0){
		//no derivatives
		return ih[0]*(1.+sq5r2+FIVETHIRDS*r2)*core;
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
                if (r2==0.) {return 0.;}
                int i =0;
		while (V[i]==0){i+=1;} //i no indexes he required dimension
                double x = (x1[i]-x2[i]);
		double l = ih[i+1];
                return ih[0]*l*x*(-FIVETHIRDS)*(1+sq5r2)*core*double(sign);
	}
	else if (S==2){
		if (P==3){
			//second derivative
                        int i =0;
			while (V[i]==0){i+=1;}
			double x = (x1[i]-x2[i]);
			double l = ih[i+1];
                        return ih[0]*FIVETHIRDS*(5.*pow(l*x,2)-l*(1.+sq5r2))*core*double(sign);
		}
		else if (P==4){
			//first derivative on two axes
                        if (r2==0.) {return 0.;}
			int i = 0;
			while (V[i]==0){i+=1;}
			int j = i+1;
			while (V[j]==0){j+=1;}
                        
                        
                        return ih[0]*ih[j+1]*(x1[j]-x2[j])*ih[i+1]*(x1[i]-x2[i])*5.*FIVETHIRDS*core*double(sign);
		}
		else{
			printf("invalid derivatives %d %d",d1,d2);
			return 0.;
		}
	}
	else if (S==3){
		if (P==4){
			//third derivative
                        if (r2==0.) {return 0.;}
			int i = 0;
			while (V[i]==0){i+=1;}
			double x = (x1[i]-x2[i]);
			double l = ih[i+1];
                        
                        return ih[0]*5*FIVETHIRDS*x*l*l*(5-5*l*pow(x,2)*oversq5r2)*core*double(sign);
		}
		else if (P==6){
			//second and first derivative
			//j is the repeated axis
                        if (r2==0.) {return 0.;}
			int i = 0;
			while (V[i]!=1){i+=1;}
			int j = 0;
			while (V[j]!=2){j+=1;}
                        double xi = (x1[i]-x2[i]);
			double li = ih[i+1];
			double xj = (x1[j]-x2[j]);
			double lj = ih[j+1];
                        
                        return ih[0]*5.*FIVETHIRDS*(li*lj*xi-5*oversq5r2*lj*lj*li*xj*xj*xi)*core*double(sign);
		}
		else if (P==8){
			//three first derivatives
                        if (r2==0.) {return 0.;}
                        int i = 0;
			while (V[i]==0){i+=1;}
			int j = i+1;
			while (V[j]==0){j+=1;}
			int k = j+1;
			while (V[k]==0){k+=1;}
			double xi = (x1[i]-x2[i]);
			double li = ih[i+1];
			double xj = (x1[j]-x2[j]);
			double lj = ih[j+1];
                        double xk = (x1[k]-x2[k]);
			double lk = ih[k+1];
                        return -ih[0]*xi*xj*xk*li*lj*lk*25.*FIVETHIRDS*oversq5r2*core*double(sign);
		}
		else{
			printf("invalid derivatives %d %d",d1,d2);
			return 0.;
		}
	}
	else if(S==4){
		if (P==16){
			//four first derivatives
                        if (r2==0.) {return 0.;}
			int i = 0;
			while (V[i]==0){i+=1;}
			int j = i+1;
			while (V[j]==0){j+=1;}
			int k = j+1;
			while (V[k]==0){k+=1;}
			int l = k+1;
			while (V[k]==0){l+=1;}
                        double xi = (x1[i]-x2[i]);
			double li = ih[i+1];
			double xj = (x1[j]-x2[j]);
			double lj = ih[j+1];
                        double xk = (x1[k]-x2[k]);
			double lk = ih[k+1];
                        double xl = (x1[l]-x2[l]);
			double ll = ih[l+1];
                        
                        return ih[0]*xi*xj*xk*xl*li*lj*lk*ll*5.*FIVETHIRDS*(sq5r2+5.*r2)*(1./pow(r2,2))*core*double(sign);
		}
		else if (P==12){
			//one second and two first
                        if (r2==0.) {return 0.;}
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
                        //k is the repeated axis
                        return ih[0]*5.*FIVETHIRDS*xi*xj*li*lj*lk*oversq5r2*(lk*pow(xk,2)*(1-sq5r2)*pow(oversq5r2,2)*5. -1.)*core*double(sign);
		}
		else if (P==9){
			//two second derivatives
                        
                        
			int i = 0;
			while (V[i]==0){i+=1;}
			int j = i+1;
			while (V[j]==0){j+=1;}
			double xi = (x1[i]-x2[i]);
			double li = ih[i+1];
			double xj = (x1[j]-x2[j]);
			double lj = ih[j+1];
                        if (r2==0.) {return ih[0]*5.*FIVETHIRDS*li*lj*core*double(sign);}
                        
                        return ih[0]*5.*FIVETHIRDS*li*lj*(1.-5.*oversq5r2*(pow(xi,2)*li+pow(xj,2)*lj)+25.*pow(oversq5r2*xi*xj,2)*li*lj*(1+oversq5r2))*core*double(sign);
		}
		else if (P==8){
			//third and first derivative
                        //j is the repeated axis
                        if (r2==0.) {return 0.;}
			int i = 0;
			while (V[i]!=1){i+=1;}
			int j = 0;
			while (V[j]!=3){j+=1;}
			double xi = (x1[i]-x2[i]);
			double li = ih[i+1];
			double xj = (x1[j]-x2[j]);
			double lj = ih[j+1];
			
                        return ih[0]*25.*FIVETHIRDS*pow(xj*lj,2)*xj*xi*li*oversq5r2*(lj*(1.-sq5r2)*(1./r2)-3.)*core*double(sign);
		}
		else if (P==5){
                    //fourth derivative
                    
                    int i = 0;
                    while (V[i]==0){i+=1;}
                    double xi = (x1[i]-x2[i]);
                    double li = ih[i+1];
                    if (r2==0.) {return ih[0]*5.*FIVETHIRDS*3.*li*li*core*double(sign);}
                    
                    return ih[0]*5.*FIVETHIRDS*(3.*pow(li,2)-30.*pow(xi*li,2)*li*oversq5r2+25.*pow(xi*li*oversq5r2,3)*xi*li*(1.+sq5r2))*core*double(sign);
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

double mat52cs(double *x1, double *x2, int d1, int d2, int D, double* ih, double* smodel){
    
    smodel[0] = ih[D+1];
    return mat52(x1,x2,d1,d2,D,ih,smodel);
}
            
int mat52csconv(double *h, int D, double* ih){
    ih[0] = pow(h[0],2);
    for (int i=1; i<D+1; i++){
	ih[i] = 1./pow(h[i],2);
        
    }
    ih[D+1] = h[D+1];
    return 0;
}
//all are searched in log space
int mat52cssearchconv(double *h, int D, double* ih){
    for (int i=0; i<D+2; i++){
        //printf("_%f",h[i]);
	ih[i] = pow(10.,h[i]);
    }
    
    return 0;
}