import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
import prediction
from collections import defaultdict
import pickle

class keydefaultdict(defaultdict):
    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError( key )
        else:
            ret = self[key] = self.default_factory(key)
            return ret

Lsupport = sp.stats.gamma.ppf(np.linspace(0,1,1002)[1:-1],4.,scale = 0.2)
Lprior = np.ones_like(Lsupport)/float(Lsupport.size)

def rfn(x):
    b = x[0]
    basevar=x[1]
    cfncoef = 60*basevar
    cfn = lambda v: cfncoef/v
    r = prediction.optatBcfoverL(b,cfn,Lsupport,Lprior,bnds=(-8,2.))
    return r

try:
    D = pickle.load(open('rcache.p','r'))[0]
except:
    D = keydefaultdict(rfn)
nb=80
nv=41
baxis = np.linspace(0.25*3600,20*3600,nb)
vaxis = np.logspace(-6,-2,nv)

rk = D[(3600.,1e-2)]
R = dict()
for k in rk.keys():
    R[k] = np.empty([nb,nv])
try:
    for i,b in enumerate(baxis):
        for j,v in enumerate(vaxis):
            print('i={} j={} : b={} v={}'.format(i,j,b,np.log10(v)))
            r = D[(b,v)]
            for k in R.keys():
                R[k][i,j] = r[k]

except KeyboardInterrupt:
    pass
pickle.dump([D],open('rcache.p','w'))
for k in ['Rmean','obsvar','c']:
    fig,ax = plt.subplots()
    CS = ax.contour(np.log10(vaxis),baxis,np.log10(R[k]),10)
    plt.clabel(CS, inline=1, fontsize=10)
    fig.savefig('figs/{}contour.png'.format(k))
    plt.close(fig)
for k in ['Esteps','Eover','B']:
    fig,ax = plt.subplots()
    CS = ax.contour(np.log10(vaxis),baxis,R[k],10)
    plt.clabel(CS, inline=1, fontsize=10)
    fig.savefig('figs/{}contour.png'.format(k))
    plt.close(fig)
for k in [1]:
    fig,ax = plt.subplots()
    CS = ax.contour(np.log10(vaxis),baxis,R['Eover']/R['B'])
    plt.clabel(CS, inline=1, fontsize=10)
    fig.savefig('figs/{}contour.png'.format('overhead'))
    plt.close(fig)
