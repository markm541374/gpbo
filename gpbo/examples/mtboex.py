import gpbo
import numpy as np
import scipy as sp
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--offset', dest='offset', action='store', default=0,type=int)

args = parser.parse_args()


rpath='results'
nreps=1
D=2

s=1e-5
lb = sp.array([-1.,-1.])
ub = sp.array([1.,1.])
tmpfile = open(os.path.join(rpath,'evallog_{}.txt'.format(args.offset)),'w')
def f(x,**ev):
    y=-np.cos(2*x[0])-np.cos(2*x[1])-ev['xa']
    print('evat {} {} {}'.format(x,ev,y))
    c=1
    tmpfile.write(','.join(['{}'.format(i) for i in x])+',{} ,{} ,{}\n'.format(ev['xa'],y,c))
    tmpfile.flush()
    return np.exp(-np.cos(2*x[0])-np.cos(2*x[1])-ev['xa']),ev['xa']+1.,{}

#from objective import truemin
all2confs=[]
all3confs=[]

#---------------
#fabolas
C={'ninit':20,
   'nsteps':21}
#all3confs.append(['fabolas',C])

#fabolas
C={'ninit':20,
   'nsteps':25,
   'switchestimator':True,
   'switchkernel':True,
   'timelimit':60*10}

all3confs.append(['fabmod2',C])


gpbo.runexp(f,lb,ub,rpath,nreps,all3confs,indexoffset=args.offset*nreps)
tmpfile.close()