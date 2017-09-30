from gpbo.core import objectives
import numpy as np
import scipy as sp
from gpbo.core import GPdc
from scipy.optimize import minimize
from matplotlib import pyplot as plt
D =2
lb = sp.array([-1.]*D)
ub = sp.array([1.]*D)
f, xm, truemin = objectives.genmat52ojf(D,lb,ub,A=1.,ls=[0.25,1.5],fixs=-1,ki=GPdc.SQUEXP)
f2, xm, truemin = objectives.genmat52ojf(D,lb,ub,A=1.,ls=[0.25,0.25],fixs=-1,ki=GPdc.SQUEXP)
g = lambda x,y: f(np.array([x,y]),**{'s':1e-9,'d':{sp.NaN}})[0]
g2 = lambda x,y: f2(np.array([x,y]),**{'s':1e-9,'d':{sp.NaN}})[0]
x =  np.linspace(-1,1,100)

Z = np.empty([100,100])
Z2 = np.empty([100,100])
for i in range(100):
    for j in range(100):
        Z[i,j]=g(x[i],x[j])
        Z2[i,j]=g2(x[i],x[j])
f,a = plt.subplots(nrows=1,ncols=2)
a[0].contour(x,x,Z)
a[1].contour(x,x,Z2)
plt.show(block=True)
