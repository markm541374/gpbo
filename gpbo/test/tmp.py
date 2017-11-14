import numpy as np
import scipy as sp
import gpbo
import pandas as pd

df = pd.read_csv('/home/mark/gpbo/exps/stopping/colville/results/s_2.csv')
dim=4
X = np.vstack([df[' x{}'.format(i)].values for i in range(dim)]).T
Y = np.array([df[' y'].values]).T
print(Y.min())
S = np.array([[0.]*Y.size]).T
D=[[np.NaN]]*Y.size
O = gpbo.core.optimize.optstate()
n=50
for i in range(n):
    O.update(X[i,:],{'d':[np.NaN],'s':0.},Y[i,0],1.,1.)
#C = gpbo.core.config.eihypdefault(lambda x:0,2,10,0.,'','')
C=gpbo.core.config.switchdefault(lambda x:0,dim,10,250,0.,'','')
persist = {'n':n,'d':dim,'flip':False,'raiseS':False}
#C.aqpara['choosereturn']={'offsetEI':0}
#gpbo.core.acquisitions.eihypaq(O,persist,**C.aqpara)
gpbo.core.choosers.globallocalregret(O,persist,**C.choosepara)
print(0)