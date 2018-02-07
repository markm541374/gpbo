import gpbo
from gpbo.core import objectives
import numpy as np
import scipy as sp
#mode='run'

gpbo.core.debugoutput=True
gpbo.core.debugoptions={'datavis':False,'drawlap':False,'cost1d':False,'ctaq':False,'support':False,'adaptive':True,'logstate':False}
mode='plot'
nreps=1
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--offset', dest='offset', action='store', default=0,type=int)

args = parser.parse_args()

vers=[2,3][0]
D=3

s=0.
lb = sp.array([-1.]*D)
ub = sp.array([1.]*D)


f = objectives.shifthart3
truemin =0
all2confs=[]
all3confs=[]
rpath='results2d'
#-----------------------

#all2confs.append(['switching_no',None])

all2confs.append(['switching_2',None])

all2confs.append(['switching_4',None])

all2confs.append(['switching_6',None])


keys ={ 'switching_2':'Switching: $R_{global}=10^{-2}$',
        'switching_4':'Switching: $R_{global}=10^{-4}$',
        'switching_6':'Switching: $R_{global}=10^{-6}$',
        'eihyp_6':'EI: $PI_{stop}=10^{-6}$',
        'eihyp_12':'EI: $PI_{stop}=10^{-12}$',
        'pesfs_1':'PES: $AQ_{stop}=10^{-1}$',
        'pesfs_4':'PES: $AQ_{stop}=10^{-4}$',
        'switching_direct':'DIRECT',
        'switching_cmaes':'CMAES'}
lfn = lambda x:keys[x]
if mode=='run':
    if vers==2:
        gpbo.runexp(f,lb,ub,rpath,nreps,all2confs,indexoffset=args.offset*nreps)
    else:
        gpbo.runexp(f,lb,ub,rpath,nreps,all3confs,indexoffset=args.offset*nreps)
elif mode=='plot':
    gpbo.plotall(all2confs+all3confs,10,rpath,trueopt=truemin+1e-99,logx=False,labelfn=lfn,showends=True,needed=[20],legend=False)
else:
    pass
