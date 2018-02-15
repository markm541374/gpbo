import numpy as np
import scipy as sp
import os
from matplotlib import pyplot as plt
import gpbo
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

path = 'results'
files = os.listdir(path)
base = 'eihyp_3_500_'

fig0,a0 = plt.subplots(1)
fig1,a1 = plt.subplots(1)
fig2,a2 = plt.subplots(1)
for i,noise in enumerate([-2,-4,-6,-8]):
    print(noise)
    names = [f for f in files if f.startswith(base+str(1000*noise))]

    D = [gpbo.optimize.readoptdata(os.path.join(path,n)) for n in names]
    [a0.plot(D[k]['index'],D[k]['trueyatxrecc'],colors[i]) for k in range(len(D))]
    [a1.plot(np.cumsum(D[k]['c']),D[k]['trueyatxrecc'],colors[i]) for k in range(len(D))]
    [a2.plot(np.cumsum(D[k]['c']+D[k]['taq']),D[k]['trueyatxrecc'],colors[i]) for k in range(len(D))]

a0.set_title('Expected Improvement')
a0.set_xlabel('Steps')
a0.set_ylabel('Immediate Regret')
a0.set_xscale('log')
a0.set_yscale('log')
a0.legend()
fig0.savefig('figs/alliter.png')

a1.set_yscale('log')
a1.set_xscale('log')
a1.set_ylabel('Immediate Regret')
a1.set_xlabel('Evaluation Cost')
#a1.legend()
#a1.set_xlim(0,1e5)
fig1.savefig('figs/allev.png')

a2.set_yscale('log')
a2.set_xscale('log')
a2.set_ylabel('Immediate Regret')
a2.set_ylabel('Total Cost')
#a2.legend()
fig2.savefig('figs/allaq.png')
