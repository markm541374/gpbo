import numpy as np
import scipy as sp
import os
from matplotlib import pyplot as plt
import gpbo
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
from gpbo import datapath
path = os.path.join(datapath,'exps/predictive/overhead/fixedPES/results')
files = os.listdir(path)
base = 'pes_3_500_'

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
a0.set_title('Predictive Entropy Search')
a0.set_xlabel('Steps')
a0.set_ylabel('Immediate Regret')
a0.set_xscale('log')
a0.set_yscale('log')
#a0.legend()
fig0.savefig('figs/iterpes.pdf')

a1.set_yscale('log')
a1.set_xscale('log')
a1.set_ylabel('Immediate Regret')
a1.set_xlabel('Evaluation Cost')
#a1.legend()
#a1.set_xlim(0,1e5)
fig1.savefig('figs/evcostpes.pdf')

a2.set_yscale('log')
a2.set_xscale('log')
a2.set_ylabel('Immediate Regret')
a2.set_xlabel('Total Cost')
#a2.legend()
fig2.savefig('figs/aqcostpes.pdf')
