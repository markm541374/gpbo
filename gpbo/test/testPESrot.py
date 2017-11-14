import numpy as np
import scipy as sp
import gpbo
import pandas as pd
gpbo.core.debugoutput['support']=True
gpbo.core.debugoutput['acqfn2d']=True
dim=2
X = np.random.uniform(-1,1,size=[30,2])
Y = np.vstack([(x[0]-0.7)**2+2*x[1]**2 for x in X])
S = np.array([[0.]*Y.size]).T
D=[[np.NaN]]*Y.size
O = gpbo.core.optimize.optstate()
for i in range(Y.size):
    O.update(X[i,:],{'d':[np.NaN],'s':0.},Y[i,0],1.,1.)
C = gpbo.core.config.pesfspredictive(lambda x:0,2,10,0.,'','')

C.aqpara['DH_SAMPLES'] =5
C.aqpara['DM_SAMPLES'] = 80
t  = sp.pi/6.
RotM = sp.array([[sp.cos(t),-sp.sin(t)],[sp.sin(t),sp.cos(t)]])
persist = {'n':30,'d':dim,'flip':False,'raiseS':False,'R':RotM}
C.aqpara['choosereturn']={'reuseH':[np.array([1.,0.5,0.5])]}
xo,_,__,___ = gpbo.core.acquisitions.PESfsaq(O,persist,**C.aqpara)
xo2,_,__,___ = gpbo.core.acquisitions.PESfsaq(O,persist,**C.aqpara)
persist['R']=np.eye(2)
xi,_,__,___ = gpbo.core.acquisitions.PESfsaq(O,persist,**C.aqpara)

print(xo,xo2,xi)
print(0)