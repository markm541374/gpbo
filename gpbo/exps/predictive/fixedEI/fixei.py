import numpy as np
import scipy as sp
import argparse
import gpbo

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--index', dest='index', action='store', default=0,type=int)
args=parser.parse_args()

kindex = gpbo.core.GPdc.MAT52
lengthscale=0.5
D=3
fpath = 'results'

lb = [-1]*D
ub = [1]*D

fc, xm, truemin = gpbo.core.objectives.genmat52ojf(D,lb,ub,A=1.,ls=lengthscale,fixs=-1,ki=kindex)

def cfn(sigma):
    return 30*(1e-6)/sigma
def icfn(c):
    return 30*(1e-6)/c
def f(x,**ev):
    y,c,aux = fc(x,**ev)
    return y,cfn(ev['s']),aux

detrandominit = np.random.uniform(-1,1,[10,3])
def detrandomaq(optstate,persist,**para):
    return detrandominit[optstate.n,:],para['ev'],persist,dict()

def aq(optstate,persist,**para):
    if optstate.n<10:
        return detrandomaq(optstate,persist,**para)
    else:
        return gpbo.core.acquisitions.eihypaq(optstate,persist,**para)

class pesfsgamma(gpbo.core.config.pesfsdefault):
    def __init__(self,*args,**kwargs):
        super(pesfsgamma,self).__init__(*args,**kwargs)
        D = len(self.aqpara['lb'])
        self.reccpara['kindex']=self.aqpara['kindex']= gpbo.core.GPdc.MAT52
        self.reccpara['mprior']=self.aqpara['mprior']= sp.array([2.]+[3.]*D)
        self.reccpara['sprior']=self.aqpara['sprior']= sp.array([0.5]+[0.15]*D)
        self.reccpara['priorshape']=self.aqpara['priorshape']='gamma'

vrange = np.linspace(-4,-8,3)
for i in range(len(vrange)):
    fname = 'eihyp_{}_{}_{}_{}.csv'.format(D,int(lengthscale*1000),int(1000*vrange[i]),args.index)
    C=pesfsgamma(f,D,10,10**vrange[i],fpath,fname)
    C.aqfn = aq
    C.stopfn = gpbo.optimize.nstopfn
    C.stoppara['nmax']=12
    C.aqpara['costfn']=cfn
    C.aqpara['icostfn']=icfn

    out = gpbo.search(C)
