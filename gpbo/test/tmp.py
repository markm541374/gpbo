import numpy as np
import scipy as sp
import gpbo
import pandas as pd

df = pd.read_csv('/home/mark/gpbo/exps/stopping/colville/results/s_1.csv')

X = np.vstack([[df[' x0'].values,df[' x1'].values]]).T
Y = np.array([df[' y'].values]).T
print(Y.min())
S = np.array([[0.]*Y.size]).T
D=[[np.NaN]]*Y.size
O = gpbo.core.optimize.optstate()
for i in range(30):
    O.update(X[i,:],{'d':[np.NaN],'s':0.},Y[i,0],1.,1.)
C = gpbo.core.config.eihypdefault(lambda x:0,2,10,0.,'','')
persist = {'n':30,'d':2}
#C.aqpara['choosereturn']={'offsetEI':0}
gpbo.core.acquisitions.eihypaq(O,persist,**C.aqpara)
print(0)