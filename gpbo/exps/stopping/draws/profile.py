import gpbo
from gpbo.core import objectives
import numpy as np
import scipy as sp
#mode='run'



all2confs=[]
rpath='friday'
#-----------------------



conf4 = [['switching_no',{'full':False,'oracle':True,'N':map(int,sp.linspace(1,100,3))}],
         ['switching_4',{'full':True,'oracle':False,'N':[]}],
         ['eihyp_nos',{'full':False,'oracle':True,'N':map(int,sp.linspace(1,100,3))}],
         ['switching_direct',{'full':False,'oracle':True,'N':map(int,sp.linspace(1,5000,3))}],
        ['switching_cmaes',{'full':False,'oracle':True,'N':map(int,sp.linspace(1,5000,3))}]]

conf6 = [['switching_no',{'full':False,'oracle':True,'N':map(int,sp.linspace(1,100,3))}],
         ['switching_6',{'full':True,'oracle':False,'N':[]}],
         ['eihyp_nos',{'full':False,'oracle':True,'N':map(int,sp.linspace(1,100,3))}],
         ['switching_direct',{'full':False,'oracle':True,'N':map(int,sp.linspace(1,5000,3))}],
        ['switching_cmaes',{'full':False,'oracle':True,'N':map(int,sp.linspace(1,5000,3))}]]

gpbo.opts.plotprofile(conf4,10,rpath,tol=0.9,target=1e-4)
gpbo.opts.plotprofile(conf6,10,rpath,tol=0.9,target=1e-6)
#import os

#D = gpbo.optimize.readoptdata(os.path.join(rpath,'{}_{}.csv'.format('switching_no',0)))
