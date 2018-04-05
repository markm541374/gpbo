import numpy as np
import scipy as sp
import pandas as pd
import sys
import time
from scipy.special import gamma
import gpbo
import gpbo.core.GPdc as GPdc
from sklearn.svm import SVR
t0 = time.clock()
#df = pd.read_csv('data/airfoilselfnoise.csv',sep='\t').as_matrix()
#df = pd.read_csv('data/yachthydro.csv',sep='\s+').as_matrix()
#df = pd.read_csv('data/Concrete_Data.csv',sep=',').as_matrix()
#df = pd.read_csv('data/ee1.csv',sep=',').as_matrix()
df = pd.read_csv('data/combinecycle.csv',sep=',').as_matrix()
t1 = time.clock()
print( 'read time {0:e}'.format(t1 - t0))
N = df.shape[0]
dim = df.shape[1]-1
print( '{}D data with {} points'.format(dim,N))

#fc, xm, truemin = gpbo.core.objectives.genmat52ojf(4,[-1]*4,[1]*4,A=1.,ls=0.9,fixs=-1,ki=GPdc.MAT52)
#X = np.random.uniform([-1],[1],size=[150,4])
#Y = np.vstack([fc(X[i,:],**{'d':[[sp.NaN]],'s':0})[0] for i in range(X.shape[0])])
#Yn=sp.random.normal(scale=0.1,size=Y.shape)
#Y+=Yn
#N = X.shape[0]
#dim=X.shape[1]
#print(Y)
np.random.seed(1234)
m = np.random.permutation(range(N))
X = df[m[:int(0.9*N)],:-1]
Y = df[m[:int(0.9*N)],-1:]
Xt = df[m[int(0.9*N):],:-1]
Yt = df[m[int(0.9*N):],-1:]
def f(x, **ev):
    print("called ojf: "+str(x))
    t0 = time.clock()
    clf = SVR(C=(4*x[0]+4+1e-6)**2, epsilon=10**(5*x[1]-5))
    clf.fit(X, Y)
    normsquerr = 1-clf.score(Xt,Yt)
    t1 = time.clock()
    print('result {}, extime {}'.format(normsquerr,t1-t0))
    return normsquerr, t1-t0,dict()

if __name__=="__main__":
    f(np.array([1.,-0.9]),**{})
    q = lambda x:f(np.array(x),**{})[0]
    sp.optimize.minimize(q,[0.99,-0.9])

    #GPdc.searchMLEhyp(X,Y,np.ones_like(Y)*1e-1,[[np.NaN]]*Y.size,np.array([-2.]*(dim+1)),np.array([2.]*(dim+1)),GPdc.SQUEXP)