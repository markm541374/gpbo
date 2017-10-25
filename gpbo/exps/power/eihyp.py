import gpbo
import numpy as np
import scipy as sp


n=100
s=0.0005

from objective import f,D

C=gpbo.core.config.pesfsdefault(f,D,100,s,'results','eihyp.csv')
C.stopfn = gpbo.core.optimize.nstopfn
#C.reccfn = gpbo.core.argminrecc
C.reccpara['check']=True
C.aqpara['nrandinit']=C.reccpara['onlyafter']=20
#print(f([-1,1,-1,1] ,**{}))
#Y = np.empty(400)
#for i in range(400):
#    Y[i] = f([0.5]*D,**{})[0]
#print(Y)
#print(np.mean(Y),np.var(Y))
out = gpbo.search(C)
print(out)

