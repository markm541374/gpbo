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
s=0.

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
        aq0 = gpbo.core.acquisitions.EIMAPaq
        aq0para = {
            'ev': {'s': s, 'd': [sp.NaN]},
            'lb': [-1.] * D,
            'ub': [1.] * D,
            'nrandinit': 10,
            'mprior': sp.array([1.] + [0.] * D),
            'sprior': sp.array([1.] * (D + 1)),
            'kindex': GPdc.MAT52,
            'volper': 1e-6
        }

        aq1 = gpbo.core.acquisitions.splocalaq
        aq1para = {
            'ev': {'s': s, 'd': [sp.NaN]},
            'lb': [-1.] * D,
            'ub': [1.] * D,
            'start':[0.]*D
        }

        choosepara = {
            'ev': aq0para['ev'],
            'lb': aq0para['lb'],
            'ub': aq0para['ub'],
            'mprior': aq0para['mprior'],
            'sprior': aq0para['sprior'],
            'kindex': aq0para['kindex'],
            'nhyp' : 12,
            'maxf': 20000,
            'onlyafter': aq0para['nrandinit'],
            'check': True,
            'everyn': 1,
            'support': 1500,
            'draws': 8000,
            'starts': 20,
            'cheatymin': ymin,
            'cheatf': f,
            'budget': n
        }
        self.aqfn = gpbo.core.acquisitions.choiceaq
        self.aqpara = {
            'aqoptions': [[aq0, aq0para], [aq1, aq1para]],
            'chooser': gpbo.core.choosers.gpcommitment,
            'choosepara': choosepara,
            'ev': aq0para['ev'],

        }


        self.stoppara = {'nmax': n}
        self.stopfn = gpbo.core.optimize.nstopfn


        self.reccfn = gpbo.core.reccomenders.gpmaprecc
        self.reccpara = {
            'ev': aq0para['ev'],
            'lb': aq0para['lb'],
            'ub': aq0para['ub'],
            'mprior': aq0para['mprior'],
            'sprior': aq0para['sprior'],
            'kindex': aq0para['kindex'],
            'volper': 1e-6,
            'onlyafter': aq0para['nrandinit'],
            'check': True,
            'everyn': 1
        }
        self.ojfchar = {'dx': len(aq0para['lb']), 'dev': len(aq0para['ev'])}
        self.ojf = f

        self.path = path
        self.fname = fname
        return

C = conf(f,D,80,s,'results','adaptive.csv')

out = gpbo.search(C)

