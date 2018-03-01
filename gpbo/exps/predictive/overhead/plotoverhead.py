import os
import dill as pickle
import numpy as np
from matplotlib import pyplot as plt

import gpbo.exps.predictive.overhead.overheads as overheads

rpath = 'plotfigs'
fpath = 'plotfigs'
fnames=[i for i in os.listdir(rpath) if i.startswith('eihyp_3_500'.format(i)) ]

T = overheads.getT(rpath, fnames, 250)
if False:
    M = overheads.buildmodel(T)
    def OM(X):
        r= overheads.cvmodel(X, M)
        return r
    pickle.dump([OM,M],open('overmodel.p','w'))
else:
    OM,M = pickle.load(open('overmodel.p','r'))
s = 60
B = 550*s

cols = plt.rcParams['axes.prop_cycle'].by_key()['color']

f0,a0 = plt.subplots()
f1,a1 = plt.subplots()
for i in np.random.randint(0,T.shape[1],size=18):
    a0.plot(T[:,i],cols[0],linewidth=0.5)
    a1.plot(np.cumsum(T[:,i]),cols[0],linewidth=0.5)
a0.plot([],[],cols[0],linewidth=0.5,label='Overhead samples')
Xp = np.arange(550)#T.shape[0])
pm,pv = OM(Xp)
std = np.sqrt(np.array([np.sum(pv[:i,:i]) for i in range(pv.shape[0])]))

a0.plot(Xp,pm,cols[1],label='Model mean')
a0.plot(Xp,pm+2*np.sqrt(np.diagonal(pv)),cols[1],linestyle='--',label=u'Model $\pm2$ standard deviation')
a0.plot(Xp,pm-2*np.sqrt(np.diagonal(pv)),cols[1],linestyle='--')


a1.plot(B-s*np.arange(int(B/s)),cols[2],linestyle=':')
a1.plot(Xp,np.cumsum(pm),cols[1])
a1.plot(Xp,np.cumsum(pm)+2*std,cols[1],linestyle='--')
a1.plot(Xp,np.cumsum(pm)-2*std,cols[1],linestyle='--')

a0.plot([],[],cols[3],label='Terminal step probability')
a0.plot([],[],cols[2],linestyle=':',label='Available overhead budget')

a1t = a1.twinx()
a1t.grid(False)
a1t.plot(overheads.stepprobs(s, B, M), cols[3])

a0.set_xlabel('Steps')
a1.set_xlabel('Steps')
a0.set_ylabel('Per-step overhead time (s)')
a1.set_ylabel('Cumulative overhead time (s)')
a1t.set_ylabel('Terminal Step Probability')
a0.legend()
f0.savefig(os.path.join(fpath,'overheadsingle.pdf'))
f1.savefig(os.path.join(fpath,'overheadcum.pdf'))
