import gpbo
from gpbo.core import objectives
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
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


keys =[ 'BLOSSOM: $R_{global}=10^{-2}$',
        'BLOSSOM: $R_{global}=10^{-4}$',
        'BLOSSOM: $R_{global}=10^{-6}$',
        'EI: $PI_{stop}=10^{-6}$',
        'EI: $PI_{stop}=10^{-12}$',
        'PES: $AQ_{stop}=10^{-1}$',
        'PES: $AQ_{stop}=10^{-4}$',
        'DIRECT',
        'CMAES']
colorlist = plt.rcParams['axes.prop_cycle'].by_key()['color']
red_patch = mpatches.Patch(color='red', label='The red data')
legendelements = [red_patch]
legendelements = [mlines.Line2D([],[],color=colorlist[i], label=k) for i,k in enumerate(keys)]
#mlines.Line2D([], [], color='blue', marker='*',
     #                     markersize=15, label='Blue stars')
lfn = lambda x:'asdf'
if mode=='run':
    if vers==2:
        gpbo.runexp(f,lb,ub,rpath,nreps,all2confs,indexoffset=args.offset*nreps)
    else:
        gpbo.runexp(f,lb,ub,rpath,nreps,all3confs,indexoffset=args.offset*nreps)
elif mode=='plot':
    gpbo.plotall(all2confs+all3confs,10,rpath,trueopt=truemin+1e-99,logx=False,labelfn=lfn,showends=True,needed=[20],legend='override',legendreplace=legendelements)
else:
    pass
