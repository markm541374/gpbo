# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.
import scipy as sp
cimport numpy as np
from scipy import linalg as spl
from scipy import stats as sps
from libc.math cimport isnan
def PhiR(x):
    # return sp.exp(sps.norm.logpdf(x) - sps.norm.logcdf(x))
    return sps.norm.pdf(x)/sps.norm.cdf(x)

def gaussian_fusion(m1,m2,V1,V2):
    V = V1.dot(spl.solve(V1+V2,V2))
    m = V.dot((spl.solve(V1,m1)+spl.solve(V2,m2)).T)
    return m,V

def expectation_prop(m0,V0,Y,Z,F,z):
    cdef int i
    needed = [True]*V0.shape[0]
    for i in xrange(V0.shape[0]):
        needed[i] =  not [Z[i]*(m0[0,i]-Z[i]*5*V0[i,i])>Z[i]*Y[i]]
    if any(needed):
        print "EP not needed for all values (>5std): "+str(needed)
    try:
        return expectation_prop_inner(m0,V0,Y,Z,F,z,needed)
    except:
        import pickle
        pickle.dump([m0,V0,Y,Z,F,z],open('epkill.p','wb'))
        raise
def expectation_prop_inner(m0,V0,Y,Z,F,z,needed):
    #expectation propagation on multivariate gaussian for soft inequality constraint
    #m0,v0 are mean vector , covariance before EP
    #Y is inequality value, Z is sign, 1 for geq, -1 for leq, F is softness variance
    #z is number of ep rounds to run
    #returns mt, Vt the value and variance for observations created by ep
    cdef int i,j,k
    m0=sp.array(m0).flatten()
    V0=sp.array(V0)
    cdef int n = V0.shape[0]
    #print "expectation prpagation running on "+str(n)+" dimensions for "+str(z)+" loops:"
    mt =sp.zeros(n)
    Vt= sp.eye(n)*float(1e10)
    m = sp.empty(n)
    V = sp.empty([n,n])
    conv = sp.empty(z)
    cdef double alpha,pr,beta,kappa,v_,m_,tmp,delta
    for i in xrange(z):
        
        #compute the m V give ep obs
        m,V = gaussian_fusion(m0,mt,V0,Vt)
        mtprev=mt.copy()
        Vtprev=Vt.copy()
        for j in [k for k in range(n) if needed[k]]:
            print [i,j]
            #the cavity dist at index j
            tmp = 1./(Vt[j,j]-V[j,j])
            v_ = (V[j,j]*Vt[j,j])*tmp
            m_ = tmp*(m[j]*Vt[j, j]-mt[j]*V[j, j])
            alpha = sp.sign(Z[j])*(m_-Y[j]) / (sp.sqrt(v_+F[j]))
            pr = PhiR(alpha)

            if isnan(pr):
                pr = -alpha
            beta = pr*(pr+alpha)/(v_+F[j])
            kappa = sp.sign(Z[j])*(pr+alpha) / (sp.sqrt(v_+F[j]))
            
            #print [alpha,beta,kappa,pr]
            mt[j] = m_+1./kappa
            #mt[j] = min(abs(mt[j]),1e5)*sp.sign(mt[j])
            Vt[j,j] = min(1e10,1./beta - v_)
        #print sp.amax(mtprev-mt)
        #print sp.amax(sp.diagonal(Vtprev)-sp.diagonal(Vt))
        #TODO make this a ratio instead of absolute
        delta = max(sp.amax(mtprev-mt),sp.amax(sp.diagonal(Vtprev)-sp.diagonal(Vt)))
        conv[i]=delta
    #print "EP finished with final max deltas "+str(conv[-3:])
    V = V0.dot(spl.solve(V0+Vt,Vt))
    m = V.dot((spl.solve(V0,m0)+spl.solve(Vt,mt)).T)
    return mt, Vt
