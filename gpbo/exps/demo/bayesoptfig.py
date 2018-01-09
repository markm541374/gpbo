import numpy as np
import scipy as sp
import gpbo
from matplotlib import pyplot as plt
f,ax = plt.subplots(nrows=4,sharex=True,figsize=[10,6])

cmap = plt.rcParams['axes.prop_cycle'].by_key()['color']
n = 200
Z = np.linspace(-1,1,n).reshape([n,1])
F = lambda x:(1+x**2)*np.sin(x*2*np.pi)
I = F(Z)
X = np.array([[-0.95,-0.25,0.2,0.6]]).T
Y = F(X)
D = [[np.NaN]]*Y.size
S = 1e-3*np.ones_like(Y)
twins=[]
for i in range(4):
    a = ax[i]
    k = gpbo.core.GPdc.kernel(gpbo.core.GPdc.MAT52,1,[1.,0.4])
    g = gpbo.core.GPdc.GPcore(X,Y,S,D,k)

    m,V = g.infer_diag_post(Z, [[np.NaN]]*n)

    def plotbounds(a,m,s):
        a.plot(Z.flatten(),m,color=cmap[0],label=u'\\mathcal{GP} Prediction')
        for i in range(3):
            a.fill_between(Z.flatten(),m-(i+1)*s,m+(i+1)*s,facecolor=cmap[0],edgecolor=None,alpha=0.2)

    plotbounds(a,m.flatten()[:n],np.sqrt(V.flatten()[:n]))

    a.plot(Z,I,'k--',linewidth=0.7,label='True Objective')

    E = g.infer_EI(Z,[[np.NaN]]*n)
    a2 = a.twinx()
    twins.append(a2)
    a2.grid(False)
    a2.plot(Z,E.T,color=cmap[3],linewidth=1.2)

    for i in range(Y.size):
        if i==0:
            a.plot(X[i],Y[i],color=cmap[1],marker='o',linestyle='none',label='Function Observations')
        else:
            a.plot(X[i],Y[i],color=cmap[1],marker='o')

    inew = np.argmax(E)
    newx = Z[inew]
    a2.plot(newx,E.flatten()[inew],'.',color=cmap[3],label='EI Acquisition Function')
    X = np.vstack([X,newx])
    Y = np.vstack([Y,F(newx)])
    D = D+[[np.NaN]]
    S = np.vstack([S,1e-3])
    a.plot([],[],color=cmap[3],label='EI Acquisition Funciton')
    a.set_ylim([-2.5,3.])
ax[0].legend(loc=9,ncol=4,bbox_to_anchor=[0.5,1.07])
ax[0].set_xlim([-1,1])
#twins[1].set_ylabel(u'EI')
#ax[1].set_ylabel(u'$y$')

f.text(0.01, 0.5, 'Function Value', ha='center', va='center', rotation='vertical')
f.text(0.99, 0.5, 'Expected Improvement', ha='center', va='center', rotation='vertical')
ax[-1].set_xlabel(u'$x$')
plt.tight_layout()
f.savefig('figs/bayesoptdemo.pdf')

