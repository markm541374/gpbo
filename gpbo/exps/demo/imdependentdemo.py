import numpy as np
import scipy as sp
import gpbo
from matplotlib import pyplot as plt
f,ax = plt.subplots(nrows=1,sharex=True,figsize=[10,3.5])
ax=[ax]
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
for i in range(1):
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
        #a2.plot(Z[m0],E.T[m0],color=cmap[3],linewidth=1.2)
        #a2.plot(Z[m1],E.T[m1],color=cmap[3],linewidth=1.2)
        #a2.plot(Z[m2],E.T[m2],linestyle='--',color=cmap[3],linewidth=1.2)

    for j in range(Y.size):
        if j==0:
            a.plot(X[j],Y[j],color=cmap[1],marker='o',linestyle='none',label='Observations')
        else:
            a.plot(X[j],Y[j],color=cmap[1],marker='o')
    #if i==0:
        #inew = -1#np.argmax(E)
        #newx = Z[inew]
        #a2.plot(newx,E.flatten()[inew],'o',color=cmap[3],label='EI Acquisition Function')
    h=0.002
    Q = np.linspace(-1,1,2/h)
    E = g.infer_EI(Q, [[np.NaN]]*Q.size).flatten()
    E = E/np.sum(E)
    E+=E*np.random.normal(loc=1,scale=0.08,size=E.size)
    E2 = E+np.maximum(0,np.random.normal(loc=0.5,size=E.size)*E.max()*0.02*np.sqrt(g.infer_diag_post(Q, [[np.NaN]]*Q.size)[1]).flatten())

    m0 = Q<-0.2
    m1 = Q>0.
    m2 = np.logical_not(np.logical_or(m0,m1))
    a2.plot(Q,E*m2/np.sum(E*m2),'r-')
    a2.plot(Q,E2*(1-m2)/np.sum(E2*(1-m2)),'g-')

    a.set_ylim([-2.,1.5])

    X = np.vstack([X,1.])
    Y = np.vstack([Y,F(1.)])
    D = D+[[np.NaN]]
    S = np.vstack([S,1e-6])

    a.plot([-0.2,-0.2],[0.2,0.5],color=cmap[1])
    a.plot([-0.,-0.],[0.2,0.5],color=cmap[1])
    a.plot([-0.2,-0.],[0.35,0.35],color=cmap[1])

ax[0].set_xlim([-1,1])
a.set_xlim([-1,1])
a2.set_xlim([-1,1])


#twins[1].set_ylabel(u'EI')
#ax[1].set_ylabel(u'$y$')

ax[-1].set_xlabel(u'$x$')
a.set_ylabel(u'$y$')
a2.set_ylabel(u'$p(x_{*})$')

plt.tight_layout()
f.savefig('figs/indepdemo.pdf')

