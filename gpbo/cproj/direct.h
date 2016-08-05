#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifndef DIRECT_
#define DIRECT_

    
struct type_intervallo {
	double		*cent, *dimen;
	double		fint, diam, maxdim, der;
	int			flagloc, flagdiv, flagcon;
	int			id;
	struct type_intervallo	*next;
};

struct type_vertice{
	struct type_intervallo	*inter;
	struct type_vertice		*next;
};

struct type_fpunt{
	double		f;
	struct type_intervallo	*punt;
};

struct {
	double		*lb, *ub, *xtemp, *xbar;
} mod_box;

struct {
	double		*vetf1, *vetf2, *xsud, *ysud;
	int			*mask;
} mod_suddividi;

typedef struct type_intervallo	intervallo;
typedef struct type_vertice		vertice;
typedef struct type_fpunt		fpunt;

// functions declarations

void direct(int n, double *lb, double *ub, int maxint, double fglob, double *xott, double *fbest, double funct(int, double*));
void deallocalistaint(intervallo *start);
void alloca_intervallo(int n, intervallo* primo);
void scalevars(int n, double *x, double *y);
void unscalevars(int n, double *y, double *x);
void ricintervallo(intervallo *root, vertice **convexhull, int *nconv);
void riduciconvexhull(vertice* convexhull, int *nelim, double eps, double toldiam, double fdir);
void suddividi(intervallo *inter, int n, int *nf, int *nint, double *xdir, double *fdir, double funct(int, double*));
void triplica(intervallo *curr, int n, int ind, double f1, double f2, 
			  int *nint, double *xdir, double *fdir);

double maxval(int n,double *dimen);
double norma(int n,double *x);

int minloc(int n, double *x, int *mask);


#endif /* DIRECT_ */