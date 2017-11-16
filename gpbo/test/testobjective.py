import gpbo
import scipy as sp
import numpy as np

from scipy.optimize import minimize
f = lambda x:gpbo.core.objectives.shifthart3(x,**{})[0]
#x0 = np.array([0.114614, 0.555649, 0.852547])*2-1
x0 = np.array([-0.77082223945184114, 0.11129778050761265, 0.70509396176578065])
print(f(x0))
res=minimize(f,x0,tol=1e-20)
print(res)
ymin = f(x0)
#for i in range(1000):
#    X = x0+np.random.normal(size=3,scale=1e-5,loc=0)
#    y=f(X)
 #   if ymin>y:
 #       print(y,X)
 ##       ymin=y
 #       Xmin=X