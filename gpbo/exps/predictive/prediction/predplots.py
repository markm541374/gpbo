import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
import prediction

import sys
print(sys.executable)
basevar=1e-4
basenum=100.

Lsupport = sp.stats.gamma.ppf(np.linspace(0,1,1002)[1:-1],4.,scale = 0.2)
Lprior = np.ones_like(Lsupport)/float(Lsupport.size)


B = np.array([1.,10,100])*3600
R=[]
fig,ax = plt.subplots(nrows=len(B),ncols=1,sharex=True,figsize=(12,9))
axt = [a.twinx() for a in ax]

for i in range(len(B)):
    b = B[i]
    cfncoef = b*basevar/basenum
    cfn = lambda v: cfncoef/v
    ax[i].set_title('$B = {} $, $c(\sigma)= \\frac{{ {:.3g} }}{{\\sigma^2}}$'.format(b,cfncoef))
    r = prediction.optatBcfoverL(b,cfn,Lsupport,Lprior,ax=ax[i],axt=axt[i],bnds=(-5,-1))
    ax[i].set_xscale('log')
    ax[i].set_yscale('log')
    R.append(r)

ax[-1].legend(loc=8,ncol=4,bbox_to_anchor=[0.5,-0.35])
fig.text(0.05, 0.5, 'Expected Regret', ha='center', va='center', rotation='vertical')
fig.text(0.95, 0.5, 'Epected Iterations', ha='center', va='center', rotation='vertical')
ax[0].set_xlim(1e-5,1e-1)
#plt.tight_layout()
fig.savefig('figs/margpredictions.pdf')
print('header\n')
for i in range(len(B)):
    s = "{:.3g} & $ \\frac{{ {:.3g} }}{{\\sigma^2}}$ & {:.3g} & {:.3g} & {:.3g} & {:.3g} \\\\".format(B[i],B[i]*basevar/basenum, R[i]['obsvar'],R[i]['Esteps'],R[i]['Rmean'],R[i]['Eover']/B[i])
    print(s)

B = np.array([1.,2.])*3600
R=[]


for i in range(len(B)):
    b = B[i]
    cfn = lambda v: B[0]*basevar/basenum/v

    r = prediction.optatBcfoverL(b,cfn,Lsupport,Lprior,bnds=(-5,-1))
    R.append(r)

print('header\n')
for i in range(len(B)):
    s = "{:.3g} & $ \\frac{{ {:.3g} }}{{\\sigma^2}}$ & {:.3g} & {:.3g} & {:.3g} & {:.3g} \\\\".format(B[i],B[0]*basevar/basenum, R[i]['obsvar'],R[i]['Esteps'],R[i]['Rmean'],R[i]['Eover']/B[i])
    print(s)
