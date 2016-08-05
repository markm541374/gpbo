# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.
import sys
sys.path.append("./..")
import OPTutils
import scipy as sp
from matplotlib import pyplot as plt
import GPdc
import ESutils
from tqdm import tqdm, tqdm_gui
runn=60
d=2
lb = sp.array([[-1.]*d])
ub = sp.array([[1.]*d])

[ojf,truexmin,trueymin] = OPTutils.gensquexpdraw(d,sp.array([-1.]*d),sp.array([1.]*d))



kindex = GPdc.SQUEXP

mprior = sp.array([0.]+[-1.]*d)
sprior = sp.array([1.]*(d+1))
maxf = 4000
volper=1e-7
s = 1e-5
ninit = 10

#para = [kindex,hlb,hub,maxf,s,ninit]
para = [kindex,mprior,sprior,volper,s,ninit]


OE = OPTutils.EIMLE(ojf,lb,ub,para)

para = dict()
para['kindex'] = GPdc.SQUEXP
para['mprior'] = sp.array([0.]+[-1.]*d)
para['sprior'] = sp.array([1.]*(d+1))
para['s'] = 1e-2
para['ninit'] = 10
#para['maxf'] = 2500
para['volper'] = 1e-7
para['DH_SAMPLES'] = 8
para['DM_SAMPLES'] = 8
para['DM_SUPPORT'] = 800
para['DM_SLICELCBPARA'] = 1.
para['SUPPORT_MODE'] = [ESutils.SUPPORT_SLICELCB,ESutils.SUPPORT_SLICEPM]
para['cfn'] = lambda x,s: ((1e-6)/s)**0.2
para['logsl'] = -8.
para['logsu'] = -2.
"""
import cProfile, pstats, StringIO
pr = cProfile.Profile()
pr.enable()
"""

OV = OPTutils.PESVS(ojf,lb,ub,para)
for i in tqdm(xrange(runn)):
    try:
        OV.step()
    except RuntimeError as e:
        print e
        break
    OE.step()

f,a = plt.subplots(7)

OE.plot(truexmin,trueymin,a,'r')
OV.plot(truexmin,trueymin,a,'g')
#OS.plot(truexmin,a,'b')



plt.show()