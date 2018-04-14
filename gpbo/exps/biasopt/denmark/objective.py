import numpy as np
import scipy as sp
import pandas as pd
import sys
import time
from scipy.special import gamma
import gpbo
import gpbo.core.GPdc as GPdc

t0 = time.clock()
df = pd.read_csv('data/denmark.csv',sep=',').as_matrix()
t1 = time.clock()
print( 'read time {0:e}'.format(t1 - t0))
N = df.shape[0]
dim = 2
print( '{}D data with {} points'.format(dim,N))

#fc, xm, truemin = gpbo.core.objectives.genmat52ojf(4,[-1]*4,[1]*4,A=1.,ls=0.9,fixs=-1,ki=GPdc.MAT52)
#X = np.random.uniform([-1],[1],size=[150,4])
#Y = np.vstack([fc(X[i,:],**{'d':[[sp.NaN]],'s':0})[0] for i in range(X.shape[0])])
#Yn=sp.random.normal(scale=0.1,size=Y.shape)
#Y+=Yn
#N = X.shape[0]
#dim=X.shape[1]
#print(Y)
X = df[:N,1:3]
Y = df[:N,-1:]
offs = sp.mean(Y)
Y -= offs
Y = Y/np.sqrt(np.var(Y))
Sbase = sp.ones([N, 1])
Dx = [[sp.NaN]] * N

#rescale to unit square
mn = np.min(X,axis=0)
mx = np.max(X,axis=0)
X = (X-mn)/(mx-mn)

pm = [4.]*(dim+1)
ps = [0.15]*(dim+1)
gammalogpdf = lambda x,shape,scale: (shape-1)*np.log(x)-x/scale-np.log(gamma(shape))-shape*np.log(scale)#-shape*np.log(scale)-np.log(gamma(shape))
priors = lambda s: gammalogpdf(s,1.,0.01)
def f(x, **ev):
    npts = max(1,int( (1-0.995*ev['xa'])*N ))
    print("called ojf at {} with xa={} ({} pts) ".format(x,ev['xa'],npts))
    hyp = [10**x[i] for i in range(dim+1)]
    t0 = time.clock()
    obsvar = 10**(3*x[dim+1]-3)
    print("hyp {}, {}".format(hyp,obsvar))
    llk = GPdc.GP_LKonly(X[:npts,:], Y[:npts,:], Sbase[:npts,:]*obsvar, Dx[:npts], GPdc.kernel(GPdc.MAT52, dim, sp.array(hyp))).plk(pm, ps,shape='gamma')
    prs = priors(obsvar)
    llk+=prs
    llknorm = llk*N/float(npts)
    llk=llknorm
    t1 = time.clock()
    if llk < -1.:
        out = sp.log(-llk) + 1.
    else:
        out = -llk

    print( "--->llk: {0} {1}    t: {2}".format(llk, out, t1 - t0))
    sys.stdout.flush()
    sys.stderr.flush()
    return out, t1-t0,dict()

if __name__=="__main__":
    f(np.array([0.1,0.1,0.1,0.1]),**{'xa':0.99})
    q = lambda x:f(np.array(x),**{})[0]
    #sp.optimize.minimize(q,[0.1,0.1,0.1,0.1,0.1,0.1])
    #sp.optimize.minimize(q,[0.1,0.1,0.1])

    #GPdc.searchMLEhyp(X,Y,np.ones_like(Y)*1e-1,[[np.NaN]]*Y.size,np.array([-2.]*(dim+1)),np.array([2.]*(dim+1)),GPdc.SQUEXP)