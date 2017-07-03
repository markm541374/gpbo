import gpbo
from gpbo.core import objectives
import numpy as np
import scipy as sp

from scipy.optimize import minimize

D = 3
lb = sp.array([-1.]*D)
ub = sp.array([1.]*D)
f, xm, truemin = objectives.genmat52ojf(D,lb,ub,ls=0.25,fixs=-1)


def wrap(x):
    return f(x,**{'d':[sp.NaN],'s':0.,'silent':True})[0]

for i in range(38):
    p = sp.random.normal(size=D)*1e-6

    res = minimize( wrap,xm+p,method='L-BFGS-B',bounds=tuple([(lb[j],ub[j]) for j in range(D)]),options={'ftol':1e-20})
    print(xm,res.x,wrap(xm),wrap(res.x))

