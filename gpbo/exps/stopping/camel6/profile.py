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
rpath='results'
#-----------------------



conf4 = [['switching_no',{'full':False,'oracle':True,'N':map(int,sp.linspace(1,100,3))}],
             ['switching_4',{'full':True,'oracle':False,'N':[]}],
             ['ei_fix',{'full':False,'oracle':True,'N':map(int,sp.linspace(1,100,3))}],
             ['switching_direct',{'full':False,'oracle':True,'N':map(int,sp.linspace(1,5000,3))}],
            ['switching_cmaes',{'full':False,'oracle':True,'N':map(int,sp.linspace(1,5000,3))}]]

conf6 = [['switching_no',{'full':False,'oracle':True,'N':map(int,sp.linspace(1,100,3))}],
         ['switching_6',{'full':True,'oracle':False,'N':[]}],
         ['ei_fix',{'full':False,'oracle':True,'N':map(int,sp.linspace(1,100,3))}],
         ['switching_direct',{'full':False,'oracle':True,'N':map(int,sp.linspace(1,5000,3))}],
        ['switching_cmaes',{'full':False,'oracle':True,'N':map(int,sp.linspace(1,5000,3))}]]

gpbo.opts.plotprofile(conf4,10,rpath,tol=0.9,target=1e-4)
gpbo.opts.plotprofile(conf6,10,rpath,tol=0.9,target=1e-6)
