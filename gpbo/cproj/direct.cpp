/*============================================================================================
!    DIRECT - C implementation of the well-known DiRect Algorithm. A 
!    Derivative-Free algorithm for bound  constrained global optimization problems 
!    proposed by Jones et al. (see Ref. below)
!    Copyright (C) 2011  G.Liuzzi, S.Lucidi, V.Piccialli
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!    Jones, D.R., Perttunen, C.D., Stuckman, B.E.: Lipschitzian optimization without 
!    the Lipschitz constant. J. Optim. Theory Appl. 79(1), 157–181 (1993)
!
==============================================================================================
*/
#include "direct.h"




void direct(int n, double *lb, double *ub, int maxint, double fglob, double *xbest, double *fbest, double funct(int, double*)){
	int i, nf, maxnf, nint, nconv, nelim;
	int halt, trovato, nmaxdimen;
	intervallo *primo, *start, *curr;
	vertice *convexhull, *currch;
	double fdir, toldiam, tolglob, eps;
	double maxdiam, mindiam, maxdimen, minder;
	double *xdir;
	char ch;

	xdir = (double *)malloc(n*sizeof(double));

	mod_box.lb    = (double *)malloc(n*sizeof(double));
	mod_box.ub    = (double *)malloc(n*sizeof(double));
	mod_box.xtemp = (double *)malloc(n*sizeof(double));
	mod_box.xbar  = (double *)malloc(n*sizeof(double));

	mod_suddividi.vetf1 = (double *)malloc(n*sizeof(double));
	mod_suddividi.vetf2 = (double *)malloc(n*sizeof(double));
	mod_suddividi.xsud  = (double *)malloc(n*sizeof(double));
	mod_suddividi.ysud  = (double *)malloc(n*sizeof(double));
	mod_suddividi.mask  = (int    *)malloc(n*sizeof(int   ));

	for(i=0;i<n;i++) {
		mod_box.lb[i]=lb[i];
		mod_box.ub[i]=ub[i];
		xbest[i] = (ub[i]+lb[i])/2.0;
		xdir[i]  = xbest[i];
		//printf("mb.lb[%d]=%f mb.ub[%d]=%f xbest[%d]=%f\n",
		//	i,mod_box.lb[i],i,mod_box.ub[i],i,xbest[i]);

	}
	printf("\n");

	primo  = (intervallo *)malloc(sizeof(intervallo));
	alloca_intervallo(n, primo);
	start = primo;
	
	/*for(i=0;i<4000000;i++){
		primo->next  = (intervallo *)malloc(sizeof(intervallo));
		primo = primo->next;
		alloca_intervallo(n, primo);
	}*/

	for(i=0;i<n;i++) {
		primo->cent[i]  = 0.5;
		primo->dimen[i] = 1.0;
	}
	primo->maxdim  = 1.0;
	primo->der     = 0.0;
	primo->diam    = norma(n,primo->dimen)/2.0;
	primo->flagloc = 0;
	primo->flagdiv = 1;
	primo->flagcon = 0;
	primo->id      = 1;

	unscalevars(n,primo->cent,xbest);
	*fbest      = funct(n,xbest);
	fdir        = *fbest;
	primo->fint = fdir;

	primo->next = NULL;

	nf          = 1;
	nint        = 1;
	nconv       = 1;
	nelim       = 1;
	toldiam     = 0.0 *sqrt((double)n)/2.0; 
	tolglob     = 1.e-4;
	eps         = 1.e-4;
	//maxnf       = 2000000000;
	halt        = 0;
	trovato     = 0;

	//printf("fbest = %f, nf = %d\n\n",*fbest,nf);

	while(!halt){
		curr     = start;
		maxdiam  = curr->diam;
		mindiam  = curr->diam;
		maxdimen = curr->dimen[0];
		for(i=1;i<n;i++) {
			if(maxdimen < curr->dimen[i]) maxdimen = curr->dimen[i];
		}
		nmaxdimen= 0;
		
		while (curr != NULL){
			if(maxdiam < curr->diam) maxdiam = curr->diam;
			if(mindiam > curr->diam) mindiam = curr->diam;
			if(maxdimen< maxval(n,curr->dimen)){
				maxdimen  = maxval(n,curr->dimen);
				nmaxdimen = 1;
			}
			else if(maxdimen == maxval(n,curr->dimen)){
				nmaxdimen = nmaxdimen + 1;
			}
			curr = curr->next;
		}
		if( ((*fbest-fglob)/(1.0 > fglob ? 1.0 : fglob)) < tolglob){
			trovato = 1;
		}
		else {
			trovato = 0;
		}
		//if( (nf > maxnf) || (nint > maxint) || trovato ){
		if( (nint > maxint) || trovato ){
			halt = 1;
			break;
		}

		convexhull = NULL;
		
		//printf("st %f\n",start->fint);

		ricintervallo(start,&convexhull,&nconv);
				
		//if(convexhull == NULL) printf("convexhull vuoto !!\n");

		riduciconvexhull(convexhull,&nelim,eps,toldiam,fdir);

		//printf("nconv = %d nelim = %d\n\n",nconv,nelim);
		//printf("cv    = %f\n\n",convexhull->inter->fint);

		currch = convexhull;
		
		for(i=1; i<=nelim; i++){
			currch = currch->next;
		}

		//printf("inizio suddivisioni\n\n");
		while (currch != NULL){
			//printf("currch id %d sta su conv.hull\n",currch->inter->id);
			if(currch->inter->flagdiv) 
			  suddividi(currch->inter,n,&nf,&nint,xdir,&fdir,funct);
			currch = currch->next;
		}
		//printf("fine   suddivisioni\n\n");
		
		unscalevars(n,xdir,xbest);
		*fbest = fdir;

	  /*!-------------------------------------------
		! dealloca la lista che memorizza gli
		! intervalli sul convex hull
		!------------------------------------------- */
		currch = convexhull;
		while (currch != NULL){
			currch->inter->flagcon = 0;
			currch = currch->next;
			convexhull->inter = NULL;
			convexhull->next  = NULL;
			free(convexhull);
			convexhull = currch;
		}
		
		//printf("fbest=%f nf=%d nconv=%d nelim=%d diam=%f DIAM=%f\n",fdir,nf,nconv,nelim,mindiam,maxdiam);
		//scanf("%c",&ch);
	}

	printf("\nDIRECT has satisfied the stopping condition.\n");
        //printf("a");
	//deallocate DIRECT structures
	deallocalistaint(start);
        //printf("b");
	free(xdir);
	free(mod_box.lb);
	free(mod_box.ub);
        //printf("c");
	free(mod_box.xtemp);
	free(mod_box.xbar);
	free(mod_suddividi.vetf1);
	free(mod_suddividi.vetf2);
        //printf("d");
	free(mod_suddividi.xsud);
	free(mod_suddividi.ysud);
	free(mod_suddividi.mask);
        //printf("e");
}

void suddividi(intervallo *curr, int n, int *nf, int *nint, double *xdir, double *fdir,	double funct(int, double*)){
	int i, j;
	int numtrue, ind1, ind2;

	numtrue = 0;
	for(i=0;i<n;i++){
		//printf("id: %d maxdim %f dimen[%d] %f\n",curr->id,curr->maxdim,i,curr->dimen[i]);
		if(curr->maxdim == curr->dimen[i]){
			for(j=0;j<n;j++) mod_suddividi.ysud[j] = curr->cent[j]; 
			mod_suddividi.ysud[i] = curr->cent[i] + 1.0*curr->dimen[i]/3.0;
			unscalevars(n,mod_suddividi.ysud,mod_suddividi.xsud);
			mod_suddividi.vetf1[i] = funct(n,mod_suddividi.xsud);

			mod_suddividi.ysud[i] = curr->cent[i] - 1.0*curr->dimen[i]/3.0;
			unscalevars(n,mod_suddividi.ysud,mod_suddividi.xsud);
			mod_suddividi.vetf2[i] = funct(n,mod_suddividi.xsud);
			mod_suddividi.mask[i] = 1;
			numtrue = numtrue + 1;

			*nf = *nf+2;
		}
		else{
			mod_suddividi.vetf1[i] = 1.e+30;
			mod_suddividi.vetf2[i] = 1.e+30;
			mod_suddividi.mask[i]  = 0;
		}
	}

	//printf("numtrue = %d\n",numtrue);

	for(i=1;i<=numtrue;i++){ 
		ind1 = minloc(n,mod_suddividi.vetf1,mod_suddividi.mask);
		ind2 = minloc(n,mod_suddividi.vetf2,mod_suddividi.mask);
		if( mod_suddividi.vetf1[ind1] < mod_suddividi.vetf2[ind2] ){
			mod_suddividi.mask[ind1] = 0;
			triplica(curr,n,ind1,mod_suddividi.vetf1[ind1],
								 mod_suddividi.vetf2[ind1],nint,xdir,fdir);
		}
		else{
			mod_suddividi.mask[ind2] = 0;
			triplica(curr,n,ind2,mod_suddividi.vetf1[ind2],
								 mod_suddividi.vetf2[ind2],nint,xdir,fdir);
		}
		*nint = *nint + 2;
	}
}

void triplica(intervallo *primo, int n, int ind, double f1, double f2, 
			  int *nint, double *xdir, double *fdir){
	int i;
	intervallo *secondo, *terzo;
	
	secondo = (intervallo *)malloc(sizeof(intervallo));
	terzo   = (intervallo *)malloc(sizeof(intervallo));

	alloca_intervallo(n, secondo);
	alloca_intervallo(n, terzo);

	terzo->next   = primo->next;
	secondo->next = terzo;
	primo->next   = secondo;

	for(i=0;i<n;i++){
		secondo->cent[i] = primo->cent[i];
		terzo->cent[i]   = primo->cent[i];
	}

	secondo->cent[ind] = secondo->cent[ind] + 1.0*primo->dimen[ind]/3.0;
	terzo->cent[ind]   = terzo->cent[ind]   - 1.0*primo->dimen[ind]/3.0;

	for(i=0;i<n;i++){
		secondo->dimen[i] = primo->dimen[i];
		terzo->dimen[i]   = primo->dimen[i];
	}

	primo->dimen[ind]   = primo->dimen[ind]/3.0;
	secondo->dimen[ind] = secondo->dimen[ind]/3.0;
	terzo->dimen[ind]   = terzo->dimen[ind]/3.0;

	primo->maxdim = maxval(n,primo->dimen);
	primo->diam   = norma(n,primo->dimen)/2.0;

	secondo->maxdim = maxval(n,secondo->dimen);
	secondo->diam   = norma(n,secondo->dimen)/2.0;

	terzo->maxdim = maxval(n,terzo->dimen);
	terzo->diam   = norma(n,terzo->dimen)/2.0;

	secondo->fint = f1;

	if(f1 < *fdir){
		*fdir = f1;
		for(i=0;i<n;i++) xdir[i] = secondo->cent[i];
	}

	terzo->fint = f2;

	if(f2 < *fdir){
		*fdir = f2;
		for(i=0;i<n;i++) xdir[i] = terzo->cent[i];
	}

	secondo->flagloc = 0;
	terzo->flagloc   = 0;

	secondo->flagdiv = 1;
	terzo->flagdiv   = 1;

	secondo->flagcon = 0;
	terzo->flagcon   = 0;

	secondo->id      = *nint+1;
	terzo->id        = *nint+2;

	//printf("secondo: %f id: %d cent: %f %f \n",f1,secondo->id,secondo->cent[0],secondo->cent[1]);
	//printf("  terzo: %f id: %d cent: %f %f \n",f2,  terzo->id,  terzo->cent[0],  terzo->cent[1]);

	secondo = NULL;
	terzo   = NULL;
}

void riduciconvexhull(vertice* convexhull, int *nelim, double eps, double toldiam, double fmin){
	vertice *currch;
	double L;
	int halt, nugua;
	
	*nelim = 0;
	nugua  = 0;
	halt   = 0;
	currch = convexhull;
	
	/*L=3.0;
	printf("abs(L)=%f",(L>0.0 ? L :-L));
	scanf("%s");*/

	//printf("rid.conv.hull: eps= %f fmin=%f\n",eps,fmin);

	while (!halt){
		
		if(currch != NULL){
			if(currch->inter->diam < toldiam){
				*nelim = *nelim + 1;
				currch = currch->next;
				continue;
			}
		}
		else{
			halt = 1;
			break;
		}
		
		if(currch->next != NULL){
			//printf("fint=%f fint=%f \n",(((currch->next)->inter)->fint), ((currch->inter)->fint) );
			//printf("diam=%f diam=%f \n",(((currch->next)->inter)->diam), ((currch->inter)->diam) );
			if((currch->next->inter->diam - currch->inter->diam) > 0.0){
				L = ( (((currch->next)->inter)->fint) - ((currch->inter)->fint) ) / 
					( (((currch->next)->inter)->diam) - ((currch->inter)->diam) );
				//printf("fint=%f diam=%f L=%f\n",(currch->inter)->fint, (currch->inter)->diam, L );
				if( currch->inter->fint - L*currch->inter->diam >  
					fmin - eps*(fmin>0.0 ? fmin :-fmin)  ){
					
					*nelim = *nelim + 1 + nugua;
					nugua  = 0;
					currch = currch->next;
				}
				else{
					halt = 1;
				}
			}
			else{
				nugua = nugua + 1;
				currch = currch->next;
			}
		}
		else{
			halt = 1;
		}
	}

	currch = NULL;
}

void ricintervallo(intervallo *root, vertice **convexhull, int *nconv){
	intervallo *curr, *primo;
	vertice *currch;
	double	maxdiam, minfunc, maxcos, coseno, norm2;
	int	halt;

	*convexhull = (vertice *)malloc(sizeof(vertice));

	curr  = root;
	primo = root;

	maxdiam = curr->diam;
	minfunc = curr->fint;

  /*!----------------------------------------
    ! search for the interval with max. diam.
	! and min. objective function
	!---------------------------------------- */
	while (curr->next != NULL){
		if( ( (curr->next)->fint < minfunc) ||
			(((curr->next)->diam > maxdiam) && ((curr->next)->fint - minfunc <= 1.e-9)) ) {
			primo   =  curr->next;	
			maxdiam =  primo->diam;
			minfunc =  primo->fint;
		}
		curr = curr->next;
	}

  /*!--------------------------------------
	! record the first interval belonging
	! to the convex hull so far identified
	!-------------------------------------- */
	(*convexhull)->inter = primo;
	/*printf("ch fint=%f diam=%f\n",((*convexhull)->inter)->fint,((*convexhull)->inter)->diam);
	printf("%f %f\n",pow(2.0,2.0),sqrt(16.0));
	printf("pr %f\n",primo->fint);
	printf("rt %f\n",root->fint);*/
	primo->flagcon = 1;
    currch = *convexhull;

	*nconv = 1;
	halt  = 0;
	
	while (!halt){
		//printf("nconv=%d\n\n",*nconv);
	  /*!-------------------------------------
		! among those intervals in the upper
		! right region, find the one with
		! maximum cosine with vector (1,0)
		!------------------------------------- */
		curr   = root;
		maxcos = -1.0;
		while (curr != NULL){
			//printf("id: %d cent: %f %f \n",curr->id,curr->cent[0],curr->cent[1]);
			if((curr->diam >= (currch->inter)->diam) && (!(curr->flagcon))){
				norm2 = sqrt(pow(curr->diam - (currch->inter)->diam,2.0) +
					         pow(curr->fint - (currch->inter)->fint,2.0)  );
				if(norm2 > 0.0){
					coseno = (curr->diam - (currch->inter)->diam) / norm2;
					if(coseno > maxcos){
						maxcos = coseno;
						primo  = curr;
					}
				}
				else{ // puo' aggiungere anche punti che coincidono sul piano f-d
					maxcos = 1.0;
					primo  = curr;
					break; 
				}
			}
			curr = curr->next;
		}
		if(maxcos > 0.0){
			//printf("aggiungo un int. al convexhull %f\n\n",maxcos);
			(currch->next) = (vertice *)malloc(sizeof(vertice));
			(currch->next)->inter = primo;
			currch = currch->next;
			*nconv = *nconv + 1;
			primo->flagcon = 1;
			//printf("ch fint=%f diam=%f\n",(currch->inter)->fint,(currch->inter)->diam);
		}
		else{
			halt = 1;
		}
	}

	currch->next = NULL;
	
	//currch = convexhull;
	//if(*convexhull != NULL) printf("convexhull NON vuoto in ricintervallo\n\n");
}

double maxval(int n,double *dimen){
	int i;
	double ris;

	ris = dimen[0];
	for(i=1;i<n;i++) {
		if(ris < dimen[i]) ris = dimen[i];
	}
	return ris;
}

void scalevars(int n, double *x, double *y){
	int i;

	for(i=0;i<n;i++) {
		y[i] = (x[i]-mod_box.lb[i])/(mod_box.ub[i]-mod_box.lb[i]);
	}
}

void unscalevars(int n, double *y, double *x){
	int i;

	for(i=0;i<n;i++) {
		x[i] = mod_box.lb[i] + y[i]*(mod_box.ub[i]-mod_box.lb[i]);
	}
	
}

double norma(int n, double *x){
	int i;
	double sum;

	sum = 0.0;
	for(i=0;i<n;i++) {
		sum = sum + x[i]*x[i];
	}

	return sqrt(sum);
}

void alloca_intervallo(int n, intervallo *primo){
	primo->cent  = (double *)malloc(n*sizeof(double));
	primo->dimen = (double *)malloc(n*sizeof(double));
	primo->next  = NULL;
}

void deallocalistaint(intervallo *start){
	intervallo	*temp;

	temp = start;
	while (temp != NULL) {
		temp = temp->next;
		free(start->cent);
		free(start->dimen);
		start->next = NULL;
		free(start);
		start = temp;
	};
}

int minloc(int n, double *x, int *mask){
	int i, ind;
	double xmin;

	ind = -1;
	xmin= 1.e+30;
	for(i=0;i<n;i++) {
		if( (x[i] < xmin) && mask[i] ){
			xmin = x[i];
			ind  = i;
		}
	}
	return ind;
}

