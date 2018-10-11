import os
import scipy as sp
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
rpath='results'
fpath='figs'

titles = {'branin':'Branin','camel3':'Camel 3','camel6':'Camel 6','hart3':'Hartmann 3','hart4':'Hartmann 4'}
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
def plotjb(ojfname):
    f,a = plt.subplots(1)
    mn=np.Inf
    mx=-np.Inf
    for i in range(16):
        fname ='pes4/'+ojfname+'_pesfs4_'+str(i)+'.csv'
        df = pd.read_csv(fname,header=None)
        for j in range(df.shape[0]):
            if df[1][j]>0:
                a.plot(j,df[1][j],'.',color=colors[0])
                if j<mn:
                    mn=j
                if j>mx:
                    mx=j

    a.set_yscale('log')
    a.plot([mn,mx],[1.76,1.76],color=colors[1])
    a.plot([mn,mx],[23.58,23.58],color=colors[2])
    a.set_ylim(0.01,1e6)
    a.set_xlabel('Iteration')
    a.set_ylabel('Jarques-Bera Score')
    a.set_title(titles[o])
    f.savefig(os.path.join(fpath,'jb_{}.pdf'.format(ojfname)))
for o in ['branin','camel3','camel6','hart3','hart4']:
    plotjb(o)
