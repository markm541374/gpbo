import numpy as np
import scipy as sp
import gpbo
from gpbo.core import GPdc
from gpbo.core import objectives
import argparse
import pickle
import time
from scipy.interpolate import interp2d
from scipy.optimize import minimize
from scipy.stats import gamma as gammadist
from scipy.stats import norm
from scipy.stats import lognorm
from scipy import stats

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--index', dest='index', action='store', default=0,type=int)
parser.add_argument('-d', '--dimension', dest='dimension', action='store', default=2,type=int)
args = parser.parse_args()

D = 2#args.dimension
lb = np.array([-1.]*D)
ub = np.array([1.]*D)
#lengthscales from 0.05 to 1.5
lengthscale = [0.2,0.2]#,0.8,0.7]
#outputscale will be normalized to 1
fc, xm, truemin = objectives.genmat52ojf(D,lb,ub,A=1.,ls=lengthscale,fixs=-1,ki=GPdc.SQUEXP)

def cfn(sigma):
    return 30*(1e-6)/sigma
def icfn(c):
    return 30*(1e-6)/c

def f(x,**ev):
    y,c,aux = fc(x,**ev)
    return y,cfn(ev['s']),aux

fpath = 'results'
from gpbo.core import debugoutput



def modeldouble(x, theta, support=[10, 30, 200]):
    # theta parameters are value,grad at supports
    p, q, r = support
    # left half
    b = theta[:4]
    A = np.array([[1., p, p**2, p**3], [0, 1, 2 * p, 3 * p**2], [1., q, q**2, q**3], [0, 1, 2 * q, 3 * q**2]])
    cl = np.linalg.solve(A, b)
    # right half
    b = theta[2:6]
    A = np.array([[1., q, q**2, q**3], [0, 1, 2 * q, 3 * q**2], [1., r, r**2, r**3], [0, 1, 2 * r, 3 * r**2]])
    cr = np.linalg.solve(A, b)

    yl = cl[0] + cl[1] * x + cl[2] * x ** 2 + cl[3] * x ** 3
    yr = cr[0] + cr[1] * x + cr[2] * x ** 2 + cr[3] * x ** 3
    y = np.select([x < p, x < q, x >= q], [theta[8], yl, yr])
    return y

def extraTdouble(X, Y, P, support=[10, 30, 200]):

    theta0 = np.array([P[0][0] * P[0][2], P[1][0], P[2][0] * P[2][2], P[3][0], P[4][0] * P[4][2], P[5][0] * P[5][2],
                       P[6][0] * P[6][2], P[7][0] * P[7][2], P[8][0]])
    if Y.size == 0:
        return theta0
    if Y.size < support[0]:
        thetaopt = theta0
        thetaopt[8] = np.mean(Y)
        return thetaopt

    def llk(theta):

        pa = 0
        pa += gammadist.logpdf(theta[0], P[0][0], loc=0, scale=P[0][2])

        pa += norm.logpdf(theta[1], loc=P[1][0], scale=P[1][1])

        pa += gammadist.logpdf(theta[2], P[2][0], loc=0., scale=P[2][2])
        pa += norm.logpdf(theta[3], loc=P[3][0], scale=P[3][1])

        pa += gammadist.logpdf(theta[4], P[4][0], loc=0., scale=P[4][2])
        pa += gammadist.logpdf(theta[5], P[5][0], loc=0., scale=P[5][2])

        pa += gammadist.logpdf(theta[6], P[6][0], loc=0., scale=P[6][2])
        pa += gammadist.logpdf(theta[7], P[7][0], loc=0., scale=P[7][2])

        pa += lognorm.logpdf(theta[8], s=P[8][0], loc=P[8][1], scale=P[8][2])

        Yh = modeldouble(X, theta)

        E = Y - Yh
        L = -pa
        L -= np.sum(stats.t.logpdf(E, loc=0., scale=theta[6], df=1. / theta[7]))

        return L

    res = minimize(llk, theta0, method='Nelder-Mead')
    thetaopt = res.x
    if Y.size <= support[1]:
        thetaopt[4] = theta0[4]
        thetaopt[5] = theta0[5]
    return thetaopt

def overlearn(X,Y,P):
    Theta = extraTdouble(X,Y,P)
    f = lambda x:modeldouble(x,Theta)
    return f


def eihyppredictiveaq(optstate,persist,**para):
    t0=time.clock()

    costfn = cfn
    icostfn = icfn
    x,ev,persist,aux = gpbo.core.acquisitions.eihypaq(optstate,persist,**para)
    if optstate.n<para['nrandinit']:
        persist['overhead']=time.clock()-t0
        ev['s'] = icostfn(para['B']/200.)
        return x,ev,persist,aux
    print('-------------------------------------------')
    ev['s'] = 10**-6

    lref,nref,M,V = pickle.load(open('/home/mark/Desktop/multidpara.p','rb'))[len(para['ub'])]
    P = pickle.load(open('/home/mark/Desktop/overheaddouble.p','rb'))[0]
    #overhead predictor
    opredict = overlearn(np.arange(optstate.n),np.array(optstate.aqtime),P)
    OIapprox = lambda x:np.sum([opredict(i) for i in range(x)])
    Oapprox = lambda x:opredict(x)
    #mean-var of regret model
    means = [interp2d(lref,nref,M[i].T) for i in range(len(M))]
    vars  = [interp2d(lref,nref,V[i].T) for i in range(len(V))]
    #hyperparamteres
    A = np.array([h[0] for h in aux['HYPdraws']])
    l = np.array([np.sqrt(h[1]*h[2]) for h in aux['HYPdraws']])

    #noiserange = np.linspace(-8,-4,200)
    nhyp = A.size#len(thetam[0])
    Bgone = sum(optstate.aqtime)+optstate.C
    B = para['B']
    Bremain = B-Bgone

    def nl2thetam(n,l):
        return np.array([p(l,n) for p in means])
    def thetam2r(thetam,m):
        return np.maximum(thetam[0]-thetam[1]*m,thetam[2]-thetam[3]*m)
    def predictR(nsteps,length,noise):
        #takes noise as sigma^2, passes on as log10
        theta = nl2thetam(np.log10(noise),length)
        R = thetam2r(theta,nsteps)
        return R

    #basic sanity check bounds
    sigmamin = icostfn(Bremain) #only take 1 evaluation
    cmax = costfn(sigmamin)
    acc=0
    i=0
    while acc+Oapprox(optstate.n+i)<Bremain:
        acc+=Oapprox(optstate.n+i)
        i+=1
    if i==0:
        #just go for a 1s evaluation since the overhead from this step will take us over budget
        ev['s'] = icostfn(1.)
        print('end on step {} due to budget end')
        return x,ev,persist,aux

    sigmamax = icostfn((Bremain-acc)/i)
    cmin = costfn(sigmamax)
    stepsmax=i
    print('maxcost {} at sigma={}, mincost {} at sigma={} for {} steps'.format(cmax,sigmamin,cmin,sigmamax,i))
    nlevels=100
    N = np.zeros([5,nlevels])
    N[0,:] = np.logspace(np.log10(sigmamin),np.log10(sigmamax),nlevels)
    N[1,:] = costfn(N[0,:])
    for i in range(nlevels):
        acc=0
        nsteps=0
        while acc<=Bremain:
            acc+=N[1,i]+Oapprox(optstate.n+nsteps)
            nsteps+=1
        N[2,i]=nsteps+optstate.n

    s = np.sum([np.log(e['s']) for e in optstate.ev])
    #use the logmean variance under projected policy as noise level to predict performance
    N[3,:] = np.exp(((N[2,:]-optstate.n)*np.log(N[0,:]) + s)/N[2,:])
    for i in range(nlevels):
        for j in range(nhyp):
            N[4,i] += (10**predictR(N[2,i],l[j],N[3,i]))*A[j]/float(nhyp)


    j = np.argmin(N[4,:])
    ev['s'] = N[0,j]
    expectedsteps = N[2,j]
    expectedcost = N[1,j]

    print('expected final regret {} after {} steps. Next ov prediction {} evcost {} meanl {} meanA {} lastover {}'.format(N[4,j],expectedsteps,Oapprox(optstate.n),expectedcost,np.mean(l),np.mean(A),optstate.aqtime[-1]))

    aux['lmean']=np.mean(l)
    aux['steppredict']=expectedsteps
    aux['overpredict']=Oapprox(optstate.n)
    aux['regretpredict']=N[4,j]

    if debugoutput['predictive']:
        from matplotlib import pyplot as plt
        f,a = plt.subplots(nrows=4,ncols=2,figsize=[12,8])
        #a[0,0].plot(np.sort(A),np.linspace(0,1,A.size))
        if not 'stepexpect' in persist.keys():
            persist['stepexpect']=np.array([[optstate.n,N[2,j]]])
        else:
            persist['stepexpect'] = np.vstack([persist['stepexpect'],[optstate.n,N[2,j]]])

        a[0,0].hist(l,20)
        a[0,1].hist(A,20)
        a[0,1].set_ylabel('A')
        a[0,0].set_ylabel('l')
        a[1,1].semilogx(N[0,:],N[2,:])
        a[1,1].plot(N[0,j],N[2,j],'ro')
        a[1,0].loglog(N[0,:],N[4,:])
        a[1,0].plot(N[0,j],N[4,j],'ro')
        a[3,1].loglog(N[0,:],N[3,:])
        a[3,1].loglog(N[0,j],N[3,j],'ro')


        a[2,0].semilogy(np.arange(optstate.n),[e['s'] for e in optstate.ev],'b.')
        a[2,0].plot([optstate.n,N[2,j]],[N[3,j],N[3,j]],'k')
        a[2,0].plot(N[2,j],N[3,j],'ko')
        a[2,0].plot(optstate.n,N[0,j],'ro')

        a[2,1].plot(persist['stepexpect'][:,0],persist['stepexpect'][:,1],'b.')

        a[3,0].plot(np.arange(200),Oapprox(np.arange(200)),'b')
        a[3,0].plot(np.arange(optstate.n),np.array(optstate.aqtime),'r.')
        plt.savefig('dbout/fig{}.png'.format(optstate.n))
        f.clf()
        plt.close(f)
        del (f)
    #persist['overhead']=time.clock()-t0
    aux['lmean']=np.mean(l)
    return x,ev,persist,aux

for i in range(1):
    fname = 'eihyp{}.csv'.format(i)
    C=gpbo.core.config.eihypdefault(f,D,200,1e-6,fpath,fname)
    C.reccpara['kindex']=C.aqpara['kindex']=GPdc.SQUEXP
    C.reccpara['mprior']=C.aqpara['mprior']= sp.array([0.]+[0.]*D)
    C.reccpara['sprior']=C.aqpara['sprior']= sp.array([0.5]*(D+1))

    debugoutput['predictive']=False
    C.aqpara['DH_SAMPLES']=200
    C.aqpara['B']=75*cfn(1e-6)
    C.aqfn = eihyppredictiveaq
    C.stopfn = gpbo.optimize.totalTorNstopfn
    C.stoppara['nmax']=200
    C.stoppara['tmax']=C.aqpara['B']
    out = gpbo.search(C)
