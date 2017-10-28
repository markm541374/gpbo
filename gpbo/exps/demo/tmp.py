import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
from gpbo.core import objectives
import gpbo
import time
import pandas as pd
ki = gpbo.core.GPdc.MAT52
lb = np.array([-1.,-1.])
ub = np.array([1.,1.])
Dim=2

mpri = np.array([1.,0.,0.])
spri = np.array([1.,1.,1.])
#mpri=np.array([2.]+[3.]*Dim)
#spri=np.array([0.5]+[0.15]*Dim)
f = 'results/pesfsd2.csv'
n=30
names = (open(f).readline().strip('\n')+''.join([',q{}'.format(i) for i in range(5)])).replace(' ','')
df = pd.read_csv(f,names=names.split(','),skiprows=1,engine='c')
X =np.array([[-0.2156916845,   0.270867522],
            [-0.5092975742,   0.6951208876],
            [0.5009155148 ,   -0.9134437223],
            [-0.5976178498 ,  0.4026885833],
            [-0.9668951227,   0.9581214121],
            [0.8401673065 ,   0.7174358881],
            [0.2371983011    ,-0.5844262466],
            [0.7264834009 ,   0.8503633595],
            [-0.2871198792,   0.5446639606],
            [0.8844036487 ,   0.8012410322],
            [-0.9831645473,   0.9578937077],
            [-0.9882631922,   0.9452475729],
            [-1,  0.6016731687],
            [-0.7712846417,   0.944783495],
            [-1 , -0.4605971737],
            [-0.9999531735 ,  0.8903341469],
            [-1 , 0.7599999069],
            [-1 , 0.828352858],
            [1  , -0.0428644065],
            [-1,  0.7803184196],
            [0.3792405857  ,  0.3902097931],
            [-0.9713127344 ,  0.786872744],
            [-0.5330624501 ,  -1],
            [-0.999994355  ,  -0.999994355],
            [-0.3785662875 ,  -0.4325527467],
            [-0.1866206194,   -0.9999926539],
            [0.147866879, 1],
            [-1 , -0.1401138949],
            [-0.9997128295 ,  0.818801615],
            [-0.9989703157 ,  0.7305794094]])

Y =np.array([[1.9626713596,
            2.0965281256,
            3.8730457732,
            1.9050492713,
            0.2187303415,
            2.2769432064,
            3.5477706862,
            2.3588292223,
            2.4797657925,
            2.7392670851,
            0.1951580325,
            0.1691622924,
            0.4126141175,
            0.6751742439,
            3.2402905708,
            0.0646689206,
            0.0146390648,
            0.0058519069,
            4.1487477826,
            0.0037315856,
            1.4752431441,
            0.0612721891,
            1.4460134892,
            3.3247053835,
            2.3265075169,
            2.24126755,
            2.2304136614,
            3.3452246307,
            0.0051448872,
            0.0478764645]]).T
S = np.array([[1e-6]*30]).T
D = [[np.NaN]]*30
n=30
X2 = np.hstack([np.array([df['x0'].values[:n]]).T,np.array([df['x1'].values[:n]]).T])
Y2 = df['y'].values[:n].reshape([X2.shape[0],1])
S2 = df['s'].values[:n].reshape([X2.shape[0],1])
X3 = np.vstack([df['x0'].values[:n],df['x1'].values[:n]]).T
#MAP = gpbo.core.GPdc.searchMAPhyp(np.copy(X), np.copy(Y), np.copy(S), D, mpri, spri, ki)
#H = gpbo.core.ESutils.drawhyp_plk(X,Y,S,D,ki,mpri,spri,28,chains=1,prior='gamma')
G = gpbo.core.PES.makeG(X,Y,S,D,ki, mpri,spri,28,prior='lognorm')
G2 = gpbo.core.PES.makeG(X3,Y2,S2,D,ki, mpri,spri,28,prior='lognorm')

hyp = sp.array([k.hyp for k in G.kf])
hmean = sp.mean(hyp, axis=0)
hstd = sp.sqrt(sp.var(hyp, axis=0))
hmin = hyp.min(axis=0)
hmax = hyp.max(axis=0)
hmed = sp.median(hyp,axis=0)
print('hyperparameters:\nmean {}\nmedian {}\nstd {}\nmin {}\nmax {}'.format(hmean,hmed,hstd,hmin,hmax))

hyp = sp.array([k.hyp for k in G2.kf])
hmean = sp.mean(hyp, axis=0)
hstd = sp.sqrt(sp.var(hyp, axis=0))
hmin = hyp.min(axis=0)
hmax = hyp.max(axis=0)
hmed = sp.median(hyp,axis=0)
print('hyperparameters:\nmean {}\nmedian {}\nstd {}\nmin {}\nmax {}'.format(hmean,hmed,hstd,hmin,hmax))
print(1)
