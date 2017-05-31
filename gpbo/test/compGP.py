import numpy as np
from gpbo.core import GPdc
from gpbo.core import flow
from matplotlib import pyplot as plt
import time
import random

t0=time.time()
N = 14
x = np.sort(np.random.rand(N,1)*2.-1,axis=0)
X = np.hstack([x])
D = [[np.NaN]]*N
S = np.ones([N,1])*0.25
#print X
y= (np.exp(-x)+np.sin(4*x))**2+np.random.randn(N,1)*0.5
dy = 2*(np.exp(-x)+np.sin(4*x))*(-np.exp(-x)+4*np.cos(4*x))+np.random.randn(N,1)*0.5

d2y = (2*(-np.exp(-x) + 4*np.cos(4*x))**2 + 2*(np.exp(-x) - 16*np.sin(4*x))* (np.exp(-x) + np.sin(4*x)))+np.random.randn(N,1)*0.5

Y = np.copy(y)
for i in range(N//2):
    j = random.randint(0,N-1)
    Y[j]=dy[j]
    D[j]=[0]
    j = random.randint(0,N-1)
    Y[j]=d2y[j]
    D[j]=[0,0]


#m.optimize()
def plot(m,a):
    xp = np.array([np.linspace(-1.1, 1.1, 200)]).T

    for i in range(X.shape[0]):
        if np.isnan(D[i][0]):
            a[0].plot(X[i,0], Y[i], 'r.')
        elif D[i]==[0]:
            a[1].plot(X[i,0], Y[i], 'r.')
        elif D[i]==[0,0]:
            a[2].plot(X[i,0], Y[i], 'r.')

    mean, var = m.infer_diag_post(xp,[[np.NaN]]*200)
    a[0].plot(xp, mean.flatten(), 'b')
    a[0].fill_between(xp[:,0], mean.flatten() - 2*np.sqrt(var.flatten()), mean.flatten() + 2*np.sqrt(var.flatten()), color='blue', alpha=0.2)
    a[0].twinx().semilogy(xp[:,0],var.flatten(),'purple')


    mean, var = m.infer_diag_post(xp,[[0]]*200)
    a[1].plot(xp, mean.flatten(), 'b')
    a[1].fill_between(xp[:,0], mean.flatten() - 2*np.sqrt(var.flatten()), mean.flatten() + 2*np.sqrt(var.flatten()), color='blue', alpha=0.2)
    a[1].twinx().semilogy(xp[:,0],var.flatten(),'purple')

    mean, var = m.infer_diag_post(xp,[[0,0]]*200)
    a[2].plot(xp, mean.flatten(), 'b')
    a[2].fill_between(xp[:,0], mean.flatten() - 2*np.sqrt(var.flatten()), mean.flatten() + 2*np.sqrt(var.flatten()), color='blue', alpha=0.2)
    a[2].twinx().semilogy(xp[:,0],var.flatten(),'purple')

    y= 1*(np.exp(-xp)+np.sin(4*xp))**2
    dy = 2*(np.exp(-xp)+np.sin(4*xp))*(-np.exp(-xp)+4*np.cos(4*xp))
    d2y = 1*(2*(-np.exp(-xp) + 4*np.cos(4*xp))**2 + 2*(np.exp(-xp) - 16*np.sin(4*xp))* (np.exp(-xp) + np.sin(4*xp)))
    a[0].plot(xp,y,'g--')
    a[1].plot(xp,dy,'g--')
    a[2].plot(xp,d2y,'g--')

#f,a = plt.subplots(nrows=3,ncols=2,figsize=(9, 6))
m = GPdc.GPcore(X,Y,S,D,GPdc.kernel(GPdc.SQUEXP,1,np.array([2.,0.4])))
mf = flow.GPcore(X,Y,S,D,flow.kernel(flow.SQUEXP,1,np.array([2.,0.4])))

#plot(m,a[:,0])
#plot(mf,a[:,1])

M=2

x = np.sort(np.random.rand(M,1)*2.-1,axis=0)
D = [[[np.NaN],[0],[0,0]][j] for j in np.random.random_integers(0,2,M)]

mean = m.infer_m_post(x,D)
fmean = mf.infer_m_post(x,D)
print('post_m: {}'.format(np.allclose(mean,fmean,rtol=1e-4) ))

mean,var = m.infer_diag_post(x,D)
fmean,fvar = mf.infer_diag_post(x,D)
print('post_diag: {}'.format(np.allclose(mean,fmean,rtol=1e-4) and np.allclose(var,fvar,rtol=1e-4)))

mean,var = m.infer_full_post(x,D)
fmean,fvar = mf.infer_full_post(x,D)

print('post_full: {}'.format(np.allclose(mean,fmean,rtol=1e-4) and np.allclose(var,fvar,rtol=1e-4)))

llk = [m.llk(),mf.llk()]
print('llk: {}'.format(np.allclose(llk[0],llk[1])))
#plt.show()
#print(V)
#print(time.time()-t0)
mf.m.optimize()
print(mf.m)