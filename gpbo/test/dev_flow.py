import GPflow as gpf
import numpy as np
from matplotlib import pyplot as plt
import time
import random

t0=time.time()
N = 4
x = np.sort(np.random.rand(N,1)*2.-1,axis=0)
X = np.hstack([x,0*np.ones(shape=(N,1))])
#print X
y= (np.exp(-x)+np.sin(4*x))**2+np.random.randn(N,1)*0.5
dy = 2*(np.exp(-x)+np.sin(4*x))*(-np.exp(-x)+4*np.cos(4*x))+np.random.randn(N,1)*0.5

d2y = (2*(-np.exp(-x) + 4*np.cos(4*x))**2 + 2*(np.exp(-x) - 16*np.sin(4*x))* (np.exp(-x) + np.sin(4*x)))+np.random.randn(N,1)*0.5

Y = np.copy(y)
for i in range(N//2):
    j = random.randint(0,N-1)
    Y[j]=dy[j]
    X[j,1]=1.
    j = random.randint(0,N-1)
    Y[j]=d2y[j]
    X[j,1]=2.

from gpbo.core.flow import dkern
k=dkern(1)
#k=kern.testk(1)
#k=gpf.kernels.RBF(1)
m = gpf.gpr.GPR(X,Y,kern=k)

m.kern.lengthscales=0.4
#m.kern.lengthscales.fixed=True
m.kern.variance=240.
m.likelihood.variance=0.5**2

#m.optimize()
def plot(m):
    xp = np.array([np.linspace(-1.1, 1.1, 200)]).T
    f,a = plt.subplots(3,figsize=(6, 6))

    mean, var = m.predict_f(np.hstack([xp,np.zeros(shape=xp.shape)]))
    for i in range(X.shape[0]):
        if X[i,1]==0:
            a[0].plot(X[i,0], Y[i], 'r.')
        elif X[i,1]==1:
            a[1].plot(X[i,0], Y[i], 'r.')
        elif X[i,1]==2:
            a[2].plot(X[i,0], Y[i], 'r.')

    a[0].plot(xp, mean, 'b')
    a[0].fill_between(xp[:,0], mean[:,0] - 2*np.sqrt(var[:,0]), mean[:,0] + 2*np.sqrt(var[:,0]), color='blue', alpha=0.2)
    a[0].twinx().semilogy(xp[:,0],var[:,0],'purple')


    mean, var = m.predict_f(np.hstack([xp,np.ones(shape=xp.shape)]))
    a[1].plot(xp, mean, 'b')
    a[1].fill_between(xp[:,0], mean[:,0] - 2*np.sqrt(var[:,0]), mean[:,0] + 2*np.sqrt(var[:,0]), color='blue', alpha=0.2)
    a[1].twinx().semilogy(xp[:,0],var[:,0],'purple')

    mean, var = m.predict_f(np.hstack([xp,2*np.ones(shape=xp.shape)]))
    a[2].plot(xp, mean, 'b')
    a[2].fill_between(xp[:,0], mean[:,0] - 2*np.sqrt(var[:,0]), mean[:,0] + 2*np.sqrt(var[:,0]), color='blue', alpha=0.2)
    a[2].twinx().semilogy(xp[:,0],var[:,0],'purple')

    y= 1*(np.exp(-xp)+np.sin(4*xp))**2
    dy = 2*(np.exp(-xp)+np.sin(4*xp))*(-np.exp(-xp)+4*np.cos(4*xp))
    d2y = 1*(2*(-np.exp(-xp) + 4*np.cos(4*xp))**2 + 2*(np.exp(-xp) - 16*np.sin(4*xp))* (np.exp(-xp) + np.sin(4*xp)))
    a[0].plot(xp,y,'g--')
    a[1].plot(xp,dy,'g--')
    a[2].plot(xp,d2y,'g--')
plot(m)
N = 4
z = np.linspace(-1,1,N)
Z = np.vstack([z,np.zeros(shape=(1,N))]).T
#print Z
#Z[2,1]=1.
M,V = m.predict_f_full_cov(Z)
plt.show()
print(V)
print(time.time()-t0)