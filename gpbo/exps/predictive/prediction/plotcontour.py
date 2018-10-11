import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
import prediction
from collections import defaultdict
import pickle
from multiprocessing import Pool
import tqdm
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
    D = dict()
nb=21
nv=21
#baxis = np.linspace(0.25*3600,40.25*3600,nb)
baxis = np.logspace(np.log10(3600/6.),np.log10(3600.*1000),nb)
vaxis = np.logspace(-7,-3,nv)

rk = rfn((3600.,1e-2))
R = dict()
for k in rk.keys():
    R[k] = np.empty([nb,nv])
existing = D.keys()
newlist = []
for i,b in enumerate(baxis):
    for j,v in enumerate(vaxis):
        #print('i={} j={} : b={} v={}'.format(i,j,b,np.log10(v)))
        if not (b,v) in existing:
            newlist.append((b,v))

#p = Pool(4)
#newvals = p.map(rfn,newlist)
#for i in range(len(newlist)):
#    D[newlist[i]] = newvals[i]

result = []
def do_work(X):
    i,x = X
    return i,rfn(x)
pool = Pool(processes=4)
for rv in tqdm.tqdm(pool.imap_unordered(do_work, zip(range(len(newlist)),newlist)), total=len(newlist)):
    result.append(rv)
newvals  = [x for _,x in sorted(result)]

for i in range(len(newlist)):
    D[newlist[i]] = newvals[i]
pickle.dump([D],open('rcache.p','w'))

for i,b in enumerate(baxis):
    for j,v in enumerate(vaxis):
        #print('i={} j={} : b={} v={}'.format(i,j,b,np.log10(v)))
        r = D[(b,v)]
        for k in R.keys():
            R[k][i,j] = r[k]
xeq = [np.linspace(-4.5,-2.85,20),np.linspace(-6.,-4.85,20)]
def xeqtoyeq(xeq,A):
    x0 = 10**xeq[0]
    yeq = A * 10**xeq
    return yeq
yeq = [xeqtoyeq(xeq[0],10000./1e-4),xeqtoyeq(xeq[1],10000./1e-6)]
def plotcontour(title,data,fname):
    fig,ax = plt.subplots()
    plt.plot([-6.75,-3.],[-0.75,3.],'--',dashes=(50,10),color='grey',linewidth=0.25)
    plt.plot([-7.,-5.],[1.,3.],'--',dashes=(50,10),color='grey',linewidth=0.25)
    plt.plot([-4.75,-3.],[-0.75,1.],'--',dashes=(50,10),color='grey',linewidth=0.25)
    plt.plot([-6.5],[-0.5],'k.')
    plt.annotate('A',(-6.65,-0.5))
    plt.plot([-3.02],[2.98],'k.')
    plt.annotate('B',(-3.25,2.85))
    CS = ax.contour(np.log10(vaxis),np.log10(baxis/3600.),data,10)
    plt.clabel(CS, inline=1, fontsize=10)
    ax.set_title(title)
    ax.set_xlabel('$\\mathrm{log}_{10}$ cost scale')
    ax.set_ylabel('$\\mathrm{log}_{10}$ Budget (hours)')

    fig.savefig('figs/{}contour.pdf'.format(fname))
    plt.close(fig)

plotcontour('$\\mathrm{log}_{10}$ Expected Regret',np.log10(R['Rmean']),'Rmean')
plotcontour('$\\mathrm{log}_{10}$ Predicted Optimum Variance',np.log10(R['obsvar']),'obsvar')
plotcontour('Expected Number of Steps',R['Esteps'],'Esteps')
plotcontour('Expected Overhead Fraction',R['Eover']/R['B'],'overhead')

