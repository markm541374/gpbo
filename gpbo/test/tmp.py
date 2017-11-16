import numpy as np
import scipy as sp
import gpbo
import pandas as pd
#gpbo.core.debugoutput['tmp']=True
df = pd.read_csv('/home/mark/gpbo/exps/stopping/colville/results/s_3.csv')
dim=4
X = np.vstack([df[' x{}'.format(i)].values for i in range(dim)]).T
Y = np.array([df[' y'].values]).T
print(Y.min())
S = np.array([[0.]*Y.size]).T
D=[[np.NaN]]*Y.size
O = gpbo.core.optimize.optstate()
n=10
for i in range(n):
    O.update(X[i,:],{'d':[np.NaN],'s':0.},Y[i,0],1.,1.)
C = gpbo.core.config.eihypdefault(lambda x:0,4,10,0.,'','')
#C=gpbo.core.config.switchdefault(lambda x:0,dim,10,250,0.,'','')

R=np.array([[-0.39400595, -0.02694207, -0.19059855, -0.89872444],
 [ 0.87676298, -0.04119442,  0.21425891, -0.42858232],
 [ 0.23588363, -0.43604394, -0.8634903,   0.09278525],
 [-0.1428457,  -0.89857822,  0.41490941,  0.00156942]])
persist = {'n':n,'d':dim,'flip':False,'raiseS':False}
#C.aqpara['choosereturn']={'offsetEI':0}
x,par,per,aux = gpbo.core.acquisitions.eihypaq(O,persist,**C.aqpara)
x1,par,per,aux = gpbo.core.acquisitions.eihypaq(O,persist,**C.aqpara)
persist['R']=R
xr,par,per,aux = gpbo.core.acquisitions.eihypaq(O,persist,**C.aqpara)
#x = gpbo.core.choosers.globallocalregret(O,persist,**C.choosepara)
print(x)
print(x1)
print(xr)