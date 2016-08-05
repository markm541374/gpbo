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
d=2
lb = sp.array([[-1.]*d])
ub = sp.array([[1.]*d])
sl=0.
su=1.
sfn = lambda x:1e-8
sls = (su-sl)*0.75
cfn = lambda x: sp.exp(-0.5*x.flatten()[0])
[ojf,truexmin,ymin] = OPTutils.gensquexpIPdraw(d,lb,ub,sl,su,sfn,sls,cfn)

print truexmin
print ojf(sp.array([sl,truexmin[0],truexmin[1]]),1e-8,[sp.NaN])
u = sp.empty(100)
v = sp.empty(100)
w = sp.empty(100)
sup = sp.linspace(-1,1,100)
spp = sp.linspace(sl,su,100)
f,a = plt.subplots(3)
for i in xrange(100):
    u[i] = ojf(sp.array([sl,sup[i],truexmin[1]]),1e-8,[sp.NaN])[0]
    v[i] = ojf(sp.array([sl,truexmin[0],sup[i]]),1e-8,[sp.NaN])[0]
    w[i] = ojf(sp.array([spp[i],truexmin[0],truexmin[1]]),1e-8,[sp.NaN])[0]
a[0].plot(sup,u)
a[1].plot(sup,v)
a[2].plot(spp,w)

"""
O = OPTutils.opt(ojf,lb,ub)
for i in xrange(30):
    O.step()



"""
kindex = GPdc.SQUEXP

mprior = sp.array([0.]+[-1.]*d)
sprior = sp.array([1.]*(d+1))
maxf = 4000
s = 1e-6
ninit = 10

#para = [kindex,hlb,hub,maxf,s,ninit]
para = [kindex,mprior,sprior,maxf,s,ninit]




para = dict()
para['kindex'] = GPdc.SQUEXPCS
para['mprior'] = sp.array([0.]+[-1.]*(d+1)+[-2.])
para['sprior'] = sp.array([1.]*(d+2)+[2.])
para['d'] = d+1
para['ninit'] = 10
para['volper'] = 1e-7
para['DH_SAMPLES'] = 8
para['DM_SAMPLES'] = 8
para['DM_SUPPORT'] = 800
para['DM_SLICELCBPARA'] = 1.
para['SUPPORT_MODE'] = [ESutils.SUPPORT_SLICELCB]
#para['cfn'] = lambda x,s: sp.exp(-0.5*x.flatten()[0])
para['sl'] = 0.
para['su'] = 1.
para['s'] = 0.
para['sfn'] = None
para['axis'] = 0
para['value'] = para['sl']
OI = OPTutils.PESIPS(ojf,lb,ub,para)
for i in tqdm(xrange(2)):
    OI.step()
    
    
print OI.X
print OI.R
print OI.Y
print OI.Ymin
print OI.Rreg

#f,a = plt.subplots(7)

OI.plot(sp.hstack([para['sl'],truexmin]),ymin)

#f.savefig("../../figs/testIPopt.pdf")

f,a = plt.subplots(1)
OI.plotcosts(a)
plt.show()
