import numpy as np
import scipy as sp
import gpbo
from gpbo.core import GPdc
from gpbo.core import objectives
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--index', dest='index', action='store', default=0,type=int)
parser.add_argument('-d', '--dimension', dest='dimension', action='store', default=2,type=int)
args = parser.parse_args()

D = args.dimension
lb = np.array([-1.]*D)
ub = np.array([1.]*D)
#lengthscales from 0.05 to 1.5
lengthscale = 0.2
maxit = 30
#outputscale will be normalized to 1
fc, xm, truemin = objectives.genmat52ojf(D,lb,ub,A=1.,ls=lengthscale,fixs=-1)

def cfn(sigma):
    return 30*(1e-6)/sigma

def f(x,**ev):
    y,c,aux = fc(x,**ev)
    return y,cfn(ev['s']),aux

fname = 'eihyp.csv'
fpath = 'results'
from gpbo.core import debugoutput

import pickle
import time
from scipy.interpolate import interp2d
def eihyppredictiveaq(optstate,persist,**para):
    t0=time.clock()

    x,ev,persist,aux = gpbo.core.acquisitions.eihypaq(optstate,persist,**para)
    if optstate.n<para['nrandinit']:
        persist['overhead']=time.clock()-t0
        return x,ev,persist,aux
    ev['s'] = 10**-6

    lref,nref,M,V,opara = pickle.load(open('paragrid.p','rb'))
    #print(opara)
    #oapprox = lambda x: (x >= 10) * (opara[0] * x ** opara[1] + opara[2]) + (x * opara[4] / 10.) * (x < 10)
    def oapprox(i):
        return opara[0][i]
    #print(np.vstack([range(100),[oapprox(j) for j in range(100)]]).T)
    means = [interp2d(lref,nref,M[i].T) for i in range(len(M))]
    vars  = [interp2d(lref,nref,V[i].T) for i in range(len(V))]
    print('-------------------------------------------')
    A = np.array([h[0] for h in aux['HYPdraws']])
    l = np.array([np.sqrt(h[1]*h[2]) for h in aux['HYPdraws']])
    sigma = -6
    thetam = [np.array([p(l[i],sigma) for i in range(l.size)]) for p in means]
    thetav = [np.array([p(l[i],sigma) for i in range(l.size)]) for p in vars]

    noiserange = np.linspace(-8,-4,200)
    nhyp = len(thetam[0])
    B = para['B'] #150*60
    def nl2thetam(n,l):
        return np.array([p(l,n) for p in means])
    def thetam2r(thetam,m):
        return np.maximum(thetam[0]-thetam[1]*m,thetam[2]-thetam[3]*m)
    nsteps = np.zeros_like(noiserange).astype(np.int)
    for i in range(nsteps.size):
        while B>oapprox(nsteps[i])+nsteps[i]*cfn(10**noiserange[i]) and nsteps[i]<199:
            nsteps[i]+=1
        nsteps[i]-=1
    regretnoiseacc = np.zeros(len(noiserange))
    overtotal = np.zeros(len(noiserange))
    regretnoiselist = []
    thetamlist=[]
    for i in range(len(A)):
        regretnoise = np.empty(len(noiserange))
        for j in range(noiserange.size):
            n = noiserange[j]
            length = l[i]
            thetam = nl2thetam(n,length)
            regretnoise[j] = thetam2r(thetam,nsteps[j])
            overtotal[j] = oapprox(nsteps[j])
        regretnoiseacc+=(10**regretnoise)*A[i]
        regretnoiselist.append((10**regretnoise)*A[i])

    j = np.argmin(regretnoiseacc)
    targetnoise = noiserange[j]
    ev['s'] = 10**(targetnoise)
    print([oapprox(i) for i in range(30)])
    print('expected final regret {} after {} steps. Current ov prediction {} evcost {}'.format(regretnoiseacc[j]/nhyp,nsteps[j],oapprox(optstate.n)-oapprox(optstate.n-1),cfn(ev['s'])))
    expectedsteps=nsteps[j]
    if debugoutput['predictive']:
        from matplotlib import pyplot as plt
        f,a = plt.subplots(nrows=5,ncols=2,figsize=[8,6])
        #a[0,0].plot(np.sort(A),np.linspace(0,1,A.size))
        a[0,0].plot(noiserange,overtotal/B)
        a[0,0].set_ylabel('overhead')

        a[0,1].plot(np.sort(l),np.linspace(0,1,l.size))
        a[0,1].set_ylabel('l')
        for nj in noiserange:
            mtheta = np.zeros_like(thetam)
            for i in range(len(A)):
                thetam = np.array(nl2thetam(nj,l[i]))
                mtheta+=thetam

            for i in range(len(thetam)):
                ix = i//2 +1
                iy = i%2
                a[ix,iy].plot(nj,mtheta[i]/float(len(A)),'b.')
            #a[ix,iy].plot(l,[tm[i] for tm in thetamlist],'b.',ms=8)
            #a[ix,iy].plot(l,thetami[i]+2*np.sqrt(thetavi[i]),'b,',ms=6)
            #a[ix,iy].plot(l,thetami[i]-2*np.sqrt(thetavi[i]),'b,',ms=6)
            a[ix,iy].set_ylabel(str(i))
        a[4,0].plot(noiserange,nsteps)
        a[4,0].text(targetnoise,expectedsteps,str(expectedsteps))
        for r in regretnoiselist:
            a[4,1].plot(noiserange,np.log10(r),'r',linewidth=1)
        a[4,1].plot(noiserange,np.log10(regretnoiseacc)-np.log10(nhyp))
        plt.savefig('dbout/fig{}.png'.format(optstate.n))
        f.clf()
        plt.close(f)
        del (f)
    #persist['overhead']=time.clock()-t0
    return x,ev,persist,aux

C=gpbo.core.config.eihypdefault(f,D,maxit,1e-6,fpath,fname)
#C=gpbo.core.config.eimledefault(f,D,maxit,1e-9,fpath,fname)
debugoutput['predictive']=False
C.aqpara['DH_SAMPLES']=50
C.aqpara['B']=75*cfn(1e-6)
C.aqfn = eihyppredictiveaq
C.aqpara['kindex']=GPdc.SQUEXP
C.stopfn = gpbo.optimize.totalTorNstopfn
C.stoppara['nmax']=200
C.stoppara['tmax']=C.aqpara['B']
out = gpbo.search(C)
