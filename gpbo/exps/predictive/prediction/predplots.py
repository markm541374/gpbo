import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
import prediction

basevar=1e-4
basenum=100.

Lsupport = sp.stats.gamma.ppf(np.linspace(0,1,1002)[1:-1],3.,scale = 0.2)
Lprior = np.ones_like(Lsupport)/float(Lsupport.size)

B = np.array([1.,10.])*3600
R=[]

fig,ax = plt.subplots(nrows=len(B),ncols=1,sharex=True)
axt = [a.twinx() for a in ax]

for i in range(len(B)):
    b = B[i]
    cfn = lambda v: b*basevar/basenum/v

    r = prediction.optatBcfoverL(b,cfn,Lsupport,Lprior,ax=ax[i],axt=axt[i],bnds=(-5,-1))
    R.append(r)

fig.savefig('figs/margpredictions.png')
print('header\n')
for i in range(len(B)):
    s = "{:.3g} & $ \\frac{{ {:.3g} }}{{\\sigma^2}}$ & {:.3g} & {:.3g} & {:.3g} & {:.3g} \\\\".format(B[i],B[i]*basevar/basenum, R[i]['obsvar'],R[i]['Esteps'],R[i]['Rmean'],R[i]['Eover']/B[i])
    print(s)