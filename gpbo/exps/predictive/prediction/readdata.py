import os
import numpy as np
import scipy as sp
import dill as pickle
from matplotlib import pyplot as plt
import gpbo
from gpbo import datapath

def readscenario(path):
    Lrange = np.linspace(0.05,1.25,25)
    Lfiles = (1000*Lrange).astype(np.int)
    Vrange = np.linspace(-6,-4,21)
    #Vrange = np.linspace(-4.,-2.,21)
    #Vrange = np.linspace(-3,-1,21)
    Vfiles = np.linspace(-6000.,-4000.,21).astype(np.int)
    print(Vfiles)
    base = 'eihyp_3'
    Lwts = sp.stats.gamma.pdf(Lrange,4.,scale=0.2)

    Rmu = np.empty_like(Vrange)
    Rva = np.empty_like(Vrange)
    TRmu = np.empty_like(Vrange)
    for i,v in enumerate(Vrange):
        norm=0
        R = []
        W = []
        TR = []
        for j,l in enumerate(Lrange):
            names = [f for f in os.listdir(path) if f.startswith('{}_{}_{}'.format(base,Lfiles[j],Vfiles[i]))]
            for n in names:
                cname = os.path.join(path,'cache',n[:-4]+'.p')
                if os.path.isfile(cname):
                    D = pickle.load(open(cname,'r'))[0]
                    #print('reload {}'.format(cname))
                else:
                    print('read {}'.format(cname))
                    D = gpbo.optimize.readoptdata(os.path.join(path,n))
                    pickle.dump([D],open(cname,'w'))
                if np.any(np.isnan(D['trueyatxrecc'])):
                    print(n)
                R.append(D['trueyatxrecc'].values[-1])
                W.append(Lwts[j])
                taq = np.sum(D['taq'].values)
                tec = np.sum(D['c'].values)
                TR.append(taq/(taq+tec))
            if len(R)==0:
                print(names,'{}_{}_{}'.format(base,Lfiles[j],Vfiles[i]))
        R = np.array(R)
        W = np.array(W)/np.sum(W)
        TR = np.array(TR)
        Rmu[i] = np.sum(R*W)
        Rva[i] = np.sum(W*(R-Rmu[i])**2)/R.size
        TRmu[i] = np.sum(TR*W)
        if np.isnan(Rva[i]):
            print(Rmu[i])
    #print(Vrange,Rmu,Rva,TRmu)
    return Vrange,Rmu, Rva, TRmu

path = os.path.join(datapath,'exps/predictive/prediction/scenarios/results_1h_v6')
#path = os.path.join(datapath,'exps/predictive/prediction/scenarios/results_1h_v2')
#path = os.path.join(datapath,'exps/predictive/prediction/scenarios/results_2h_v6')
Vrange,Rmu, Rva, TRmu = readscenario(path)
fig,ax = plt.subplots()
ax.plot(Vrange,Rmu)
ax.fill_between(Vrange,Rmu-2*np.sqrt(Rva),Rmu+2*np.sqrt(Rva),alpha=0.2)
ax.set_yscale('log')
ax.twinx().plot(Vrange,TRmu, 'r')
fig.savefig('figs/tmp.png')

pickle.dump([10**Vrange,Rmu,TRmu],open(os.path.join(path,'cache/out.p'),'w'))