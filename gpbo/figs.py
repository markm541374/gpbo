from __future__ import print_function
xrange=range

import scipy as sp
import numpy as np
import pandas as pd
import gpbo
import os
import matplotlib
from matplotlib import pyplot as plt
#plt.style.use('seaborn-paper')
#plt.rc('font',serif='Times')

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
def stoppingplots(path,names,n,legendnames=None,fname='',title='',offset=0.,fpath=None,showlegend=False,logy=True,r2path=False):
    if fpath is None:
        fpath=path
    if legendnames==None:
        legendnames=names
    D = dict()
    for name in names:
        D[name]=[]
        for i in range(n):
            D[name].append(gpbo.optimize.readoptdata(os.path.join(path,'{}_{}.csv'.format(name,i))))
    if r2path:
        Dno = dict()
        for name in names:
            if name.startswith('switch'):
                Dno[name]=[]
                for i in range(n):
                    d0 = D[name][i]
                    d1 = gpbo.optimize.readoptdata(os.path.join(r2path,'{}_{}.csv'.format(name,i)))
                    dc = pd.concat([d0[:int(d1['n'][0])-1],d1])
                    Dno[name].append(dc)

    f,a = plt.subplots(1)
    for j,name in enumerate(names):
        gpbo.opts.plotquartsends(a,[D[name][k]['index'] for k in range(n)],[D[name][k]['trueyatxrecc']-offset-min(0,D[name][k]['trueyatxrecc'].values[-1]) for k in range(n)],colors[j],0,legendnames[j])
    if r2path:
        for j,name in enumerate(names):
            if name.startswith('switch'):
                gpbo.opts.plotquartsends(a,[Dno[name][k]['n']-1 for k in range(n)],[Dno[name][k]['trueyatxrecc']-offset-min(0,Dno[name][k]['trueyatxrecc'].values[-1]) for k in range(n)],colors[j],0,legendnames[j],median=True,noends=True,linestyle=':',notwin=True)
    if logy:
        a.set_yscale('log')
    #a.set_ylim(10.01,10.02)
    a.set_xlabel('Steps')
    a.set_ylabel('Regret')
    a.set_title(title)
    if showlegend:
        a.legend()
    f.savefig(os.path.join(fpath,'stopping_{}.pdf'.format(fname)))

    print('tablerow:')
    E=dict()
    for j,name in enumerate(names):
        E[name]=dict()
        E[name]['steps'] = np.empty(n)
        E[name]['endRegret'] = np.empty(n)
        for i in range(n):
            E[name]['steps'][i] = D[name][i]['index'].values[-1]
            E[name]['endRegret'][i] = D[name][i]['trueyatxrecc'].values[-1]-offset-min(0,D[name][i]['trueyatxrecc'].values[-1])
        E[name]['r'] = '{:.3g}'.format(np.mean(E[name]['endRegret']))
        print(name)
        print(E[name]['endRegret'])
       # print(E[name]['endRegret']*E[name]['steps'])
        E[name]['s'] = '{:.3g}'.format(np.mean(E[name]['steps']))
        E[name]['rs']= '{:.3g}'.format(np.mean(E[name]['steps']*E[name]['endRegret']))

    print(' & '.join([k for k in names]))
    print(' & '.join([E[k]['r'] for k in names]))
    print(' & '.join([E[k]['s'] for k in names]))
    print(' & '.join([E[k]['rs'] for k in names]))

    return

def overhead(path,name0,names,n,legendnames=None,fname='',title='',fpath=None):
    if fpath is None:
        fpath=path
    if legendnames==None:
        legendnames=[name0,names]
    D = dict()
    D[name0]=[]
    for i in range(n):
        D[name0].append(gpbo.optimize.readoptdata(os.path.join(path,'{}_{}.csv'.format(name0,i))))
    for name in names:
        D[name]=[]
        for i in range(n):
            D[name].append(gpbo.optimize.readoptdata(os.path.join(path,'{}_{}.csv'.format(name,i))))
    f,a = plt.subplots(1)
    gpbo.opts.plotquarts(a,[D[name0][k]['index'] for k in range(n)],[D[name0][k]['taq'] for k in range(n)],colors[0],'-',legendnames[0])
    gpbo.opts.plotquarts(a,[D[names[0]][k]['index'] for k in range(n)],[D[names[0]][k]['taq'] for k in range(n)],colors[1],'-',legendnames[1][0])
    #for i in range(n):
    #    a.plot(D[names[0]][i]['index'],D[names[0]][i]['taq'],color=colors[1],linestyle='-')
    a.set_xlabel('Steps')
    a.set_ylabel('Overhead Time (s)')
    a.set_title(title)
    a.legend()
    f.savefig(os.path.join(fpath,'overhead_{}.pdf'.format(fname)))

    j = max([max([len(D[name][i]['taq']) for i in range(n)]) for name in names])
    print(j)
    ref = np.zeros(j)
    a =0
    for i in range(n):
        print(len(D[name0][i]['taq']))
        if len(D[name0][i]['taq']) >= j:
            ref+=D[name0][i]['taq'].values[:j]
            a+=1
    ref/=float(a)

    for name in names:

        acc=[0.,0.]
        for i in range(n):
            l = len(D[name][i]['taq'])
            tref = np.sum(ref[:l])
            taq = np.sum(D[name][i]['taq'].values[:l])
            #print(l,taq,tref)
            acc[0]+= taq
            acc[1]+=tref
        print(acc[0]/acc[1])
    return