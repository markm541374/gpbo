import numpy as np
import scipy as sp
import gpbo
import time
import tqdm
from matplotlib import pyplot as plt
ki = gpbo.core.GPdc.MAT52
Dim=3
lb = np.array([-1.]*Dim)
ub = np.array([1.]*Dim)

mpri = np.array([1.]+[0.]*Dim)
spri = np.array([1.]*(Dim+1))
def timing(n,N):
    X = np.random.uniform(-1,1,[n,Dim])
    Y = sp.vstack([gpbo.core.objectives.shifthart3(X[i,:],**{})[0] for i in range(X.shape[0])])
    S = np.ones_like(Y)
    D = [[sp.NaN]]*Y.size
    G = gpbo.core.PES.makeG(X,Y,S,D,ki, mpri,spri,24,prior='lognorm')


    Z = np.random.uniform(-1,1,[N,Dim])

    tsingle = np.zeros([N,15])
    tbatch = np.zeros([N,15])
    t0 = time.time()
    nrange = np.linspace(1,N,N,dtype=int)
    for i in tqdm.tqdm(range(N)):
        for k in range(15):
            t0 = time.time()
            a = G.infer_EI_post(Z[:(i+1),:],[[sp.NaN]]*(i+1))
            t1 = time.time()
            tbatch[i,k] = t1-t0

            t0 = time.time()
            for j in range(i+1):
                a = G.infer_EI_post(Z[j,:],[[sp.NaN]]*(j+1))
            t1 = time.time()
            tsingle[i,k] = t1-t0
    return nrange,np.median(tbatch,axis=1),np.median(tsingle,axis=1)

f,a = plt.subplots(1)
x,tb,ts  = timing(100,40)
a.plot(x,tb,color='b',linestyle='-',marker='D',label='n=100, Batch')
a.plot(x,ts,color='b',linestyle='-',marker='o',label='n=100, Sequential')

#x,tb,ts  = timing(50,30)
#a.plot(x,tb,color='g',marker='D')
#a.plot(x,ts,color='g',marker='o')
#f.savefig('figs/timing.png')

x,tb,ts  = timing(20,40)
a.plot(x,tb,color='r',marker='D',linestyle='-', label='n=20,   Batch')
a.plot(x,ts,color='r',marker='o',linestyle='-', label='n=20,   Sequential')
a.set_yscale('log')
a.set_ylabel('CPU time (s)')
a.set_xlabel('Evaluation Points')
a.legend()
f.savefig('figs/timing.png')
print(1)
