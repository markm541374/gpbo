# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

import OPTutils
import scipy as sp
from matplotlib import pyplot as plt
import GPdc
import ESutils
from tqdm import tqdm_gui
d=2
lb = sp.array([[-1.]*d])
ub = sp.array([[1.]*d])
[ojf,truexmin] = [lambda x,s,d:[-1.,-1.],None]

para = dict()
para['kindex'] = GPdc.SQUEXP
para['mprior'] = sp.array([0.]+[-1.]*d)
para['sprior'] = sp.array([1.]*(d+1))
para['s'] = 1e-6
para['ninit'] = 10
#para['maxf'] = 2500
para['volper'] = 1e-6
para['DH_SAMPLES'] = 8
para['DM_SAMPLES'] = 8
para['DM_SUPPORT'] = 400
para['DM_SLICELCBPARA'] = 0.025
para['SUPPORT_MODE'] = ESutils.SUPPORT_SLICEPM
OP = OPTutils.PESFS(ojf,lb,ub,para)


import pickle
[OP.X,OP.Y,OP.S,OP.D,OP.R,OP.C,OP.T,OP.Tr,OP.Ymin] = pickle.load(open('state.p','rb'))

OP.pcs()
#OP.step()
plt.figure(1)
#plt.plot(OP.R[-1,0],OP.R[-1,1],'r.')
plt.plot(OP.X[:,0],OP.X[:,1],'b.')

plt.show()

