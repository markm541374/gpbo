import gpbo
import numpy as np
import scipy as sp
mode='run'
#mode='plot'

D=2

s=1e-3
lb = sp.array([-1.,-1.])
ub = sp.array([1.,1.])

def f(x, **ev):
    sn=ev['xa']
    y = 3.-np.cos(2.9*(1-0.2*sn)*(x[0]-sn*0.5))-np.cos(3.1*(x[1]+sn*0.2))+sn**2
    c=0.5+3.*(1.-sn)**2

    print 'f inputs x:{} ev:{} outputs y:{}  c:{}'.format(x, ev, y ,c)
    return y , c, dict()


allconfs=[]

#-----------------------
#eimle
C=gpbo.core.config.eimledefault(f,D,12,s,'results','null.csv')
C.aqpara['nrandinit']=10
C.stoppara = {'tmax': 60*10}
C.stopfn = gpbo.core.optimize.totaltstopfn

#allconfs.append(['eimle',C])

#----------------------
#pesfs
C=gpbo.core.config.pesfsdefault(f,D,12,s,'results','null.csv')
C.stoppara = {'tmax': 60*10}
C.aqpara['nrandinit']=10
C.stopfn = gpbo.core.optimize.totaltstopfn

#allconfs.append(['pesfs',C])

#-----------------
#pesbs
C=gpbo.core.config.pesbsdefault(f,D,50,s,'results','null.csv')
C.stoppara = {'tmax': 60 * 15}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.aqpara['overhead']='predict'
C.aqpara['nrandinit']=20

#allconfs.append(['pesbs',C])
#-----------------
#mtbo
C={'lowtask':4,
   'ninit':10,
   'nsteps':50,
   'switchestimator':True}

#allconfs.append(['mtbo4_post',C])

#-----------------
#mtbo
C={'lowtask':16,
   'ninit':10,
   'nsteps':50,
   'switchestimator':True}

#allconfs.append(['mtbo16_post',C])

#-----------------
#mtbo
C={'lowtask':64,
   'ninit':10,
   'nsteps':50,
   'switchestimator':True}

#allconfs.append(['mtbo64_post',C])

#-----------------
#mtbo
C={'lowtask':4,
   'ninit':10,
   'nsteps':50,
   'switchestimator':False}

allconfs.append(['mtbo4_arg',C])
#-----------------
#mtbo
C={'lowtask':16,
   'ninit':10,
   'nsteps':50,
   'switchestimator':False}

#allconfs.append(['mtbo16_arg',C])

#-----------------
#mtbo
C={'lowtask':64,
   'ninit':10,
   'nsteps':12,
   'switchestimator':False}

#allconfs.append(['mtbo64_arg',C])

#fabolasmod------------------------------
C={'ninit':20,
   'nsteps':22,
   'switchkernel':True,
   'switchestimator':True}
allconfs.append(['fabmod',C])
#---------------
#fabolas
C={'ninit':20,
   'nsteps':22}
#allconfs.append(['fabolas',C])


if mode=='run':
    gpbo.runexp(f,lb,ub,'resultsX',1,allconfs)
elif mode=='plot':
    gpbo.plotall(allconfs,2,'results',trueopt=1.)
else:
    pass
