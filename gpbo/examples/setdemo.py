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
C.stoppara = {'tmax': 60*5}
C.stopfn = gpbo.core.optimize.totaltstopfn

#allconfs.append(['eimle',C])

#----------------------
#pesfs
C=gpbo.core.config.pesfsdefault(f,D,12,s,'results','null.csv')
C.stoppara = {'tmax': 60*1}
C.aqpara['nrandinit']=10
C.stopfn = gpbo.core.optimize.totaltstopfn

#allconfs.append(['pesfs',C])

#-----------------
#pesbs
C=gpbo.core.config.pesbsdefault(f,D,50,s,'results','null.csv')
C.stoppara = {'tmax': 60 * 5}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.aqpara['overhead']='last'
C.aqpara['nrandinit']=10

#allconfs.append(['pesbs',C])

#-----------------
#mtbo
C={'lowtask':2,
   'ninit':10,
   'nsteps':20}

#allconfs.append(['mtbo2',C])

#-----------------
#mtbo
C={'lowtask':4,
   'ninit':10,
   'nsteps':11}

#allconfs.append(['mtbo4',C])

#-----------------
#mtbo
C={'lowtask':8,
   'ninit':10,
   'nsteps':20}

#allconfs.append(['mtbo8',C])
#---------------
#fabolas
C={'ninit':12,
   'nsteps':40}
allconfs.append(['fabolas',C])
if mode=='run':
    gpbo.runexp(f,lb,ub,'results',2,allconfs)
elif mode=='plot':
    gpbo.plotall(allconfs,2,'results',trueopt=1.)
else:
    pass
