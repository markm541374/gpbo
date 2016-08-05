# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

import OPTutils
import scipy as sp
from matplotlib import pyplot as plt
import GPdc
import ESutils
from tqdm import tqdm, tqdm_gui
import copy
runn = 10
d=2
lb = sp.array([[-1.]*d])
ub = sp.array([[1.]*d])
[ojf,truexmin,trueymin] = OPTutils.gensquexpdraw(d,sp.array([-1.]*d),sp.array([1.]*d),ignores=1e-9)
print "truexmin "+ str(truexmin)
f,a = plt.subplots(8)
#ojf = OPTutils.banana
#truexmin=OPTutils.bananamin
#trueymin=0.
O = OPTutils.opt(ojf,lb,ub)
for i in tqdm(xrange(10)):
    O.step()
initstate = copy.deepcopy([O.X,O.Y,O.S,O.D,O.R,O.C,O.T,O.Tr,O.Ymin,O.Xmin,O.Yreg,O.Rreg])

    
kindex = GPdc.MAT52
mprior = sp.array([0.]+[-1.]*d)
sprior = sp.array([1.]*(d+1))
volper=1e-8
s = 1e-9
ninit = 10
para = [kindex,mprior,sprior,volper,s,ninit]

OE = OPTutils.EIMLE(ojf,lb,ub,para,initstate=copy.deepcopy(initstate))
#OL = OPTutils.LCBMLE(ojf,lb,ub,para)
#[OL.X,OL.Y,OL.S,OL.D,OL.R,OL.C,OL.T,OL.Tr,OL.Ymin] = initstate
#[OE.X,OE.Y,OE.S,OE.D,OE.R,OE.C,OE.T,OE.Tr,OE.Ymin,OE.Xmin,OE.Yreg,OE.Rreg] = copy.deepcopy(initstate)




para = dict()
para['kindex'] = GPdc.MAT52
para['mprior'] = sp.array([0.]+[-1.]*d)
para['sprior'] = sp.array([1.]*(d+1))
para['s'] = 1e-9
para['ninit'] = 10
#para['maxf'] = 2500
para['volper'] = 1e-6
para['DH_SAMPLES'] = 8
para['DM_SAMPLES'] = 8
para['DM_SUPPORT'] = 1200
para['DM_SLICELCBPARA'] = 1.
para['SUPPORT_MODE'] = [ESutils.SUPPORT_SLICELCB,ESutils.SUPPORT_SLICEPM]

OP = OPTutils.PESFS(ojf,lb,ub,para,initstate=copy.deepcopy(initstate))
#[OP.X,OP.Y,OP.S,OP.D,OP.R,OP.C,OP.T,OP.Tr,OP.Ymin,OP.Xmin,OP.Yreg,OP.Rreg] = copy.deepcopy(initstate)
for i in tqdm_gui(xrange(runn),gui=True):
    
    state = [OP.X,OP.Y,OP.S,OP.D,OP.R,OP.C,OP.T,OP.Tr,OP.Ymin]
    try:
        pass
        #OP.step()
    except:
        import pickle
        pickle.dump(state,open('state.p','wb'))
        raise
        
    OE.step()
    O.step()
    #OL.step()
    
    plt.close(f)
    f,a = plt.subplots(8)
    
    OE.plot(truexmin,trueymin,a,'r')
    OP.plot(truexmin,trueymin,a,'g')
    O.plot(truexmin,trueymin,a,'b')
    plt.draw()
    

f,a = plt.subplots(8)
    
OE.plot(truexmin,trueymin,a,'r')
OP.plot(truexmin,trueymin,a,'g')
O.plot(truexmin,trueymin,a,'b')
    #plt.draw()        


plt.show()