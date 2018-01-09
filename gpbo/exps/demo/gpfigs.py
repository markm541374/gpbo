import numpy as np
import scipy as sp
import gpbo
from matplotlib import pyplot as plt
f,ax = plt.subplots(nrows=3,ncols=2,sharex=True,figsize=[12,6])

cmap = plt.rcParams['axes.prop_cycle'].by_key()['color']
n = 200
Z = np.linspace(-1,1,n).reshape([n,1])
X = np.array([[-0.8,-0.3,0.2,0.9]]).T
Y = np.array([[-0.6,0.4,-2.6,0.5]]).T
D = [[np.NaN],[np.NaN],[0],[np.NaN]]
S = np.array([[1e-3,1e-3,1e-3,1e-3]]).T
for i in range(2):
    a = ax[:,i]
    k = gpbo.core.GPdc.kernel([gpbo.core.GPdc.SQUEXP,gpbo.core.GPdc.MAT52][i],1,[1.,0.4])
    g = gpbo.core.GPdc.GPcore(X,Y,S,D,k)

    m,V = g.infer_diag_post(np.vstack([Z]*3) , [[np.NaN]]*n + [[0]]*n + [[0,0]]*n)
    #a[0].plot(Z,m[:,:n].T,color=cmap[0])
    #a[1].plot(Z,m[:,n:2*n].T,color=cmap[0])
    #a[2].plot(Z,m[:,2*n:3*n].T,color=cmap[0])

    def plotbounds(a,m,s):
        a.plot(Z.flatten(),m,color=cmap[0])
        for i in range(3):
            a.fill_between(Z.flatten(),m-(i+1)*s,m+(i+1)*s,facecolor=cmap[0],edgecolor=None,alpha=0.2)

    plotbounds(a[0],m.flatten()[:n],np.sqrt(V.flatten()[:n]))
    plotbounds(a[1],m.flatten()[n:2*n],np.sqrt(V.flatten()[n:2*n]))
    plotbounds(a[2],m.flatten()[2*n:3*n],np.sqrt(V.flatten()[2*n:3*n]))


    m,V = g.infer_full(np.vstack([Z]*3) , [[np.NaN]]*n + [[0]]*n + [[0,0]]*n)
    W  = sp.random.multivariate_normal(m.flatten(),V,5)
    a[0].plot(Z,W[:,:n].T,color=cmap[3],linewidth='0.5')
    a[1].plot(Z,W[:,n:2*n].T,color=cmap[3],linewidth='0.5')
    a[2].plot(Z,W[:,2*n:].T,color=cmap[3],linewidth='0.5')
    #a[1].plot(Z[1:],np.diff(W.T,axis=0)*n/2.,color=cmap[3],linewidth='0.5')
    #a[2].plot(Z[1:-1],np.diff(np.diff(W.T,axis=0),axis=0)*(n/2.)**2,color=cmap[3],linewidth='0.5')

    for i in range(Y.size):
        if D[i]==[np.NaN]:
            a[0].plot(X[i],Y[i],color=cmap[1],marker='o')
        if D[i]==[0]:
            a[1].plot(X[i],Y[i],color=cmap[1],marker='o')
        if D[i]==[0,0]:
            a[2].plot(X[i],Y[i],color=cmap[1],marker='o')

ax[0,0].set_ylabel(u'$y$')
ax[1,0].set_ylabel(u'$\\frac{\partial y}{\partial x}$')
ax[2,0].set_ylabel(u'$\\frac{\partial^2 y}{\partial x^2}$')
ax[2,0].set_xlabel(u'$x$')
ax[2,1].set_xlabel(u'$x$')
ax[2,0].set_xlim([-1,1])
ax[2,1].set_xlim([-1,1])
ax[0,0].set_title(u'Squared Exponential Kernel')
ax[0,1].set_title(u'Mat\\\'ern $\\frac{5}{2}$ Kernel')
plt.tight_layout()
ax[0,0].set_ylim(ax[0,1].get_ylim())
ax[1,0].set_ylim(ax[1,1].get_ylim())
ax[2,0].set_ylim(ax[2,1].get_ylim())
f.savefig('figs/gpdemo.pdf')

