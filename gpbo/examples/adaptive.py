import matplotlib
matplotlib.use('Agg')

import gpbo
from gpbo.core import objectives as objectives
from gpbo.core import GPdc as GPdc
import os
import scipy as sp
import time

gpbo.core.debugoutput=True
gpbo.core.debugoptions={'datavis':True,'drawlap':False,'cost1d':False,'ctaq':False,'support':False,'adaptive':True}

D=2
n=50
s=1e-6

lb=[-1.,-1.]
ub = [1.,1.]
f, xmin, ymin = objectives.genmat52ojf(D,lb,ub,ls=0.2)

with open('results/matdrawopt.txt','w') as o:
    o.write('reported truemin x {} ; y {}'.format(xmin,ymin))

gpbo.core.ESutils.plot2dFtofile(f,os.path.join('dbout', 'truegeneratedobjective' + time.strftime('%d_%m_%y_%H:%M:%S') + '.png'),xmin=xmin)

class conf():
    """
        fixed s, space is [-1,1]^D

        """

    def __init__(self, f, D, n, s, path, fname):
        self.aqfn = gpbo.core.acquisitions.EIMAPaq
        self.aqpara = {
            'ev': {'s': s, 'd': [sp.NaN]},
            'lb': [-1.] * D,
            'ub': [1.] * D,
            'nrandinit': 12,
            'mprior': sp.array([1.] + [0.] * D),
            'sprior': sp.array([1.] * (D + 1)),
            'kindex': GPdc.MAT52,
            'volper': 1e-6
        }
        """
        self.aqfn = gpbo.core.acquisitions.PESfsaq
        self.aqpara = {
            'ev': {'s': s, 'd': [sp.NaN]},
            'lb': [-1.] * D,
            'ub': [1.] * D,
            'nrandinit': 12,
            'volper': 1e-6,
            'mprior': sp.array([1.] + [0.] * D),
            'sprior': sp.array([1.] * (D + 1)),
            'kindex': GPdc.MAT52,
            'DH_SAMPLES': 16,
            'DM_SAMPLES': 16,
            'DM_SUPPORT': 2000,
            'SUPPORT_MODE': [gpbo.core.ESutils.SUPPORT_LAPAPROT],
            'DM_SLICELCBPARA': 16.,
            'noS': False,
        }
        """
        self.stoppara = {'nmax': n}
        self.stopfn = gpbo.core.optimize.nstopfn

        self.reccfn = gpbo.core.reccomenders.adaptiverecc
        self.reccpara = {
            'ev': self.aqpara['ev'],
            'lb': self.aqpara['lb'],
            'ub': self.aqpara['ub'],
            'mprior': self.aqpara['mprior'],
            'sprior': self.aqpara['sprior'],
            'kindex': self.aqpara['kindex'],
            'volper': 1e-6,
            'onlyafter': self.aqpara['nrandinit'],
            'check': True,
            'everyn': 1,
            'support':5000,
            'draws':40000,
            'starts':20,
            'cheatymin':ymin,
            'cheatf':f
        }
        self.ojfchar = {'dx': len(self.aqpara['lb']), 'dev': len(self.aqpara['ev'])}
        self.ojf = f

        self.path = path
        self.fname = fname
        return

C = conf(f,D,100,s,'results','adaptive.csv')

out = gpbo.search(C)

