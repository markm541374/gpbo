import numpy as np
import scipy as sp
import os
from matplotlib import pyplot as plt
import gpbo
import dill as pickle
import gpbo
import scipy as sp

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

_,PP = pickle.load(open('/home/mark/gpbo/exps/predictive/performance/model.p'))
from gpbo.exps.predictive.performance.model import perfmodel
def PM(X,S,L):
    Xa = np.ones([X.size,2])*np.log10(S)
    Xa[:,0] = X
    truemean,log10std = perfmodel(Xa,L,PP)
    natlogmu = np.log(truemean)
    natlogvar = (np.log(10)*log10std)**2
    return natlogmu,natlogvar
from gpbo import datapath

path = os.path.join(datapath,'exps/predictive/performance/data/results16')
files = os.listdir(path)
for il,ls in enumerate([300,800,1500]):
    base = 'eihyp_3_{}_'.format(ls)

    #fig0,a0 = plt.subplots(1)
    fig1,a1 = plt.subplots(1)
    for i,noise in enumerate([-2,-4,-6,-8]):
        #print(noise)
        names = [f for f in files if f.startswith(base+str(1000*noise))]
        names = [names[k] for k in [1,3,5,6]]
        D = [gpbo.optimize.readoptdata(os.path.join(path,n)) for n in names]
        #[a0.plot(D[k]['index'][10:],D[k]['trueyatxrecc'][10:],colors[i],linewidth=0.5) for k in range(len(D))]
        [a1.plot(np.cumsum(D[k]['c'])[10:],D[k]['trueyatxrecc'].values[10:],colors[i],linewidth=0.5) for k in range(len(D))]
        xp = np.linspace(D[0]['c'].values[0]*10,D[0]['c'].values[0]*len(D[0]['c']),150)
        #low0, med0, upp0 = gpbo.core.ESutils.quartsirregular([np.cumsum(D[k]['c']) for k in range(len(D))],[D[k]['trueyatxrecc'] for k in range(len(D))],xp)
        #me,va = gpbo.core.ESutils.mvirregular([np.cumsum(D[k]['c']) for k in range(len(D))],[D[k]['trueyatxrecc'] for k in range(len(D))],xp)
        #a1.plot(xp,med0,colors[i],linewidth=0.75)
        X = np.arange(10,2*D[k]['index'].values[-1])
        mu,var = PM(X,10**noise,ls/1000.)
        #a0.plot(X,np.exp(mu),colors[i])
        #a0.fill_between(X,np.exp(mu-2*np.sqrt(var)),np.exp(mu+2*np.sqrt(var)),facecolor=colors[i],alpha=0.2)
        #a1.plot(X*D[k]['c'].values[0],np.exp(mu),colors[i],alpha=0.5)
        a1.fill_between(X*D[k]['c'].values[0],np.exp(mu-0.05*np.sqrt(var)),np.exp(mu+0.05*np.sqrt(var)),facecolor=colors[i],alpha=0.5)
        a1.fill_between(X*D[k]['c'].values[0],np.exp(mu-2*np.sqrt(var)),np.exp(mu+2*np.sqrt(var)),facecolor=colors[i],alpha=0.15)
        a1.fill_between(X*D[k]['c'].values[0],np.exp(mu-np.sqrt(var)),np.exp(mu+np.sqrt(var)),facecolor=colors[i],alpha=0.15)
        #gpbo.opts.plotquarts(a0,[D[k]['index'] for k in range(len(D))],[D[k]['trueyatxrecc'] for k in range(len(D))],colors[i],'-','$\sigma^2={}$'.format(noise))
        #gpbo.opts.plotquarts(a1,[np.cumsum(D[k]['c']) for k in range(len(D))],[D[k]['trueyatxrecc'] for k in range(len(D))],colors[i],'-',str(noise))
        #gpbo.opts.plotquarts(a2,[np.cumsum(D[k]['c']+D[k]['taq']) for k in range(len(D))],[D[k]['trueyatxrecc'] for k in range(len(D))],colors[i],'-',str(noise))
        if il==0:
            a1.plot([],[],colors[i],label="$\\sigma^2 = 10^{{{}}}$".format(noise))
    if il==0:
        a1.legend()


    #a0.set_title('Expected Improvement ls = {}'.format(ls/1000.))
    #a0.set_xlabel('Steps')
    #a0.set_ylabel('Immediate Regret')
    #a0.set_xscale('log')
    #a0.set_yscale('log')
    #a0.legend()
    #fig0.savefig('figs/iterei_{}.pdf'.format(ls))
    a1.set_xlim(2*1e-2,1e6)
    a1.set_ylim(1e-10,1e2)
    a1.set_title('L = {}'.format(ls/1000.))
    a1.set_yscale('log')
    a1.set_xscale('log')
    a1.set_ylabel('Immediate Regret')
    a1.set_xlabel('Evaluation Cost')
    fig1.savefig('figs/evcostei_{}.pdf'.format(ls))

    #a2.set_title('Expected Improvement ls = {}'.format(ls/1000.))
    #a2.set_yscale('log')
    #a2.set_xscale('log')
    #a2.set_ylabel('Immediate Regret')
    #a2.set_ylabel('Total Cost')
    #a2.legend()
    #fig2.savefig('figs/aqcostei_{}.pdf'.format(ls))

    #plt.close(fig0)
    plt.close(fig1)
    #plt.close(fig2)
