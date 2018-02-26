import numpy as np
import scipy as sp
import os
from matplotlib import pyplot as plt
import gpbo
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

path = 'data/results8'
files = os.listdir(path)
for ls in [100,200,300,500,700,900,1200]:
    base = 'eihyp_3_{}_'.format(ls)

    fig0,a0 = plt.subplots(1)
    fig1,a1 = plt.subplots(1)
    fig2,a2 = plt.subplots(1)
    for i,noise in enumerate([-2,-3,-4,-5,-6,-7,-8]):
        print(noise)
        names = [f for f in files if f.startswith(base+str(1000*noise))]

        D = [gpbo.optimize.readoptdata(os.path.join(path,n)) for n in names]
        gpbo.opts.plotquarts(a0,[D[k]['index'] for k in range(len(D))],[D[k]['trueyatxrecc'] for k in range(len(D))],colors[i],'-','$\sigma^2={}$'.format(noise))
        gpbo.opts.plotquarts(a1,[np.cumsum(D[k]['c']) for k in range(len(D))],[D[k]['trueyatxrecc'] for k in range(len(D))],colors[i],'-',str(noise))
        gpbo.opts.plotquarts(a2,[np.cumsum(D[k]['c']+D[k]['taq']) for k in range(len(D))],[D[k]['trueyatxrecc'] for k in range(len(D))],colors[i],'-',str(noise))

    X = np.logspace(-2,6,100)
    Y = 0.002*X**(-0.5)
    a1.plot(X,Y,'k--')

    X = np.logspace(1.5,2,100)
    Y = 1e12*X**(-8)
    a1.plot(X,Y,'k--')

    a0.set_title('Expected Improvement ls = {}'.format(ls/1000.))
    a0.set_xlabel('Steps')
    a0.set_ylabel('Immediate Regret')
    a0.set_xscale('log')
    a0.set_yscale('log')
    a0.legend()
    fig0.savefig('figs/iterei_{}.pdf'.format(ls))

    a1.set_title('Expected Improvement ls = {}'.format(ls/1000.))
    a1.set_yscale('log')
    a1.set_xscale('log')
    a1.set_ylabel('Immediate Regret')
    a1.set_xlabel('Evaluation Cost')
    fig1.savefig('figs/evcostei_{}.pdf'.format(ls))

    a2.set_title('Expected Improvement ls = {}'.format(ls/1000.))
    a2.set_yscale('log')
    a2.set_xscale('log')
    a2.set_ylabel('Immediate Regret')
    a2.set_ylabel('Total Cost')
    #a2.legend()
    fig2.savefig('figs/aqcostei_{}.pdf'.format(ls))

    plt.close(fig0)
    plt.close(fig1)
    plt.close(fig2)
