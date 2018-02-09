import numpy as np
import scipy as sp
import gpbo
from matplotlib import pyplot as plt
f,ax = plt.subplots(nrows=2,sharex=True,figsize=[10,7])

cmap = plt.rcParams['axes.prop_cycle'].by_key()['color']
n = 600
Z = np.linspace(-1,1,n).reshape([n,1])

#I = F(Z)
X = np.array([[-0.95,-0.75,-0.6,-0.25,-0.1,-0.18,0.,0.1,0.2,0.5,0.8]]).T
F = lambda x: -np.cos((x+0.15)*1.4*np.pi)
Y = F(X)
#Y = np.array([[1.5,2.,0.1,-0.35,0.1,2.15,1.95]]).T
D = [[np.NaN]]*Y.size
S = 1e-6*np.ones_like(Y)
twins=[]
for i in range(2):
    a = ax[i]
    k = gpbo.core.GPdc.kernel(gpbo.core.GPdc.MAT52,1,[1.,0.4])
    g = gpbo.core.GPdc.GPcore(X,Y,S,D,k)

    m,V = g.infer_diag_post(Z, [[np.NaN]]*n)

    def plotbounds(a,m,s):
        a.plot(Z.flatten(),m,color=cmap[0],label=u'\\mathcal{GP} Prediction')
        for i in range(3):
            a.fill_between(Z.flatten(),m-(i+1)*s,m+(i+1)*s,facecolor=cmap[0],edgecolor=None,alpha=0.2)

    plotbounds(a,m.flatten()[:n],np.sqrt(V.flatten()[:n]))

    #a.plot(Z,I,'k--',linewidth=0.7,label='True Objective')

    E = g.infer_EI(Z,[[np.NaN]]*n)
    if i==0:
        a2 = a.twinx()
        twins.append(a2)
        a2.grid(False)
        m0 = Z<-0.2
        m1 = Z>0.
        m2 = np.logical_not(np.logical_or(m0,m1))
        a2.plot(Z[m0],E.T[m0],color=cmap[3],linewidth=1.2)
        a2.plot(Z[m1],E.T[m1],color=cmap[3],linewidth=1.2)
        a2.plot(Z[m2],E.T[m2],linestyle='--',color=cmap[3],linewidth=1.2)

    for j in range(Y.size):
        if j==0:
            a.plot(X[j],Y[j],color=cmap[1],marker='o',linestyle='none',label='Observations')
        else:
            a.plot(X[j],Y[j],color=cmap[1],marker='o')
    if i==0:
        inew = -1#np.argmax(E)
        newx = Z[inew]
        a2.plot(newx,E.flatten()[inew],'o',color=cmap[3],label='EI Acquisition Function')
        a.plot([],[],color=cmap[3],label='EI Acquisition Funciton')
    a.set_ylim([-2.,1.5])

    X = np.vstack([X,1.])
    Y = np.vstack([Y,F(1.)])
    D = D+[[np.NaN]]
    S = np.vstack([S,1e-6])

    a.plot([-0.2,-0.2],[0.2,0.5],color=cmap[1])
    a.plot([-0.,-0.],[0.2,0.5],color=cmap[1])
    a.plot([-0.2,-0.],[0.35,0.35],color=cmap[1])
    if i==1:
        a.plot([-0.15],[-1],marker='o',linestyle='none', color=cmap[3],label='Next Evaluation')
        a.plot([],[],color=cmap[3],label='GRR Acquisition')
        a.plot([],[],color=cmap[1],label='Convex Region')

ax[1].legend(loc=8,ncol=3,bbox_to_anchor=[0.5,-0.01])
ax[0].set_xlim([-1,1])
#twins[1].set_ylabel(u'EI')
#ax[1].set_ylabel(u'$y$')

f.text(0.013, 0.5, 'Function Value', ha='center', va='center', rotation='vertical')
f.text(0.99, 0.75, 'GRR acquisition', ha='center', va='center', rotation='vertical')
ax[-1].set_xlabel(u'$x$')
plt.tight_layout()
f.savefig('figs/stoppingdemo.pdf')

