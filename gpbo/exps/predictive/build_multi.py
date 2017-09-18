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
from scipy.special import gamma

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


gammalogpdf = lambda x,shape,scale: (shape-1)*np.log(x)-x/scale#-shape*np.log(scale)-np.log(gamma(shape))
normlogpdf = lambda x,loc,scale: -((x-loc)/scale)**2
lognormlogpdf = lambda x,s,loc,scale: -((np.log(x)-loc)/scale)**2 - np.log(x)
tlogpdf = lambda x,v,scale: -0.5*(v+1)*np.log(1+((x/scale)**2)/v)

def modeldouble(x, theta, support=[10, 30, 200]):
    # theta parameters are value,grad at supports
    p, q, r = support
    # left half
    A = np.array([[1., p, p**2, p**3], [0, 1, 2 * p, 3 * p**2], [1., q, q**2, q**3], [0, 1, 2 * q, 3 * q**2]])
    cl = np.linalg.solve(A, theta[:4])
    # right half
    b = theta[2:6]
    A = np.array([[1., q, q**2, q**3], [0, 1, 2 * q, 3 * q**2], [1., r, r**2, r**3], [0, 1, 2 * r, 3 * r**2]])
    cr = np.linalg.solve(A, theta[2:6])

    yl = cl[0] + cl[1] * x + cl[2] * x ** 2 + cl[3] * x ** 3
    yr = cr[0] + cr[1] * x + cr[2] * x ** 2 + cr[3] * x ** 3
    y = theta[8]*(x < p) +yl*np.logical_and(x>=p,x<q)+yr*(x>=q)
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
        #pa += gammadist.logpdf(theta[0], P[0][0], loc=0, scale=P[0][2])
        pa += gammalogpdf(theta[0], P[0][0], P[0][2])

        #pa += norm.logpdf(theta[1], loc=P[1][0], scale=P[1][1])
        pa += normlogpdf(theta[1], P[1][0], P[1][1])

        #pa += gammadist.logpdf(theta[2], P[2][0], loc=0., scale=P[2][2])
        pa += gammalogpdf(theta[2], P[2][0], P[2][2])
        #pa += norm.logpdf(theta[3], loc=P[3][0], scale=P[3][1])
        pa += normlogpdf(theta[3], P[3][0], P[3][1])

        #pa += gammadist.logpdf(theta[4], P[4][0], loc=0., scale=P[4][2])
        pa += gammalogpdf(theta[4], P[4][0], P[4][2])
        #pa += gammadist.logpdf(theta[5], P[5][0], loc=0., scale=P[5][2])
        pa += gammalogpdf(theta[5], P[5][0], P[5][2])

        #pa += gammadist.logpdf(theta[6], P[6][0], loc=0., scale=P[6][2])
        pa += gammalogpdf(theta[6], P[6][0], P[6][2])
        #pa += gammadist.logpdf(theta[7], P[7][0], loc=0., scale=P[7][2])
        pa += gammalogpdf(theta[7], P[7][0], P[7][2])

        pa += lognormlogpdf(theta[8], P[8][0], P[8][1], P[8][2])

        Yh = modeldouble(X, theta)

        E = Y - Yh
        L = -pa
        L -= np.sum(tlogpdf(E, 1./theta[7],theta[6]))

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

def getVopt(optstate,persist,para,ev,aux):
    P = pickle.load(open('/home/mark/Desktop/overheaddouble.p','rb'))[0]
    #overhead predictor
    opredict = overlearn(np.arange(optstate.n),np.array(optstate.aqtime),P)
    OIapprox = lambda x:np.sum([opredict(i) for i in range(x)])
    Oapprox = lambda x:opredict(x)
    #hyperparamteres
    A = np.array([h[0] for h in aux['HYPdraws']])
    l_ = np.array([np.sqrt(h[1]*h[2]) for h in aux['HYPdraws']])
    ind = np.argsort(l_)
    l=l_[ind]
    A=A[ind]

    means = para['means']
    vars = para['vars']
    #noiserange = np.linspace(-8,-4,200)
    nhyp = A.size#len(thetam[0])
    Bgone = sum(optstate.aqtime)+optstate.C
    B = para['B']
    Bremain = B-Bgone

    def nl2thetam(n,l):
        return np.array([p(l,n) for p in means])

    def thetam2r(theta,x):
        phi = [10**theta[0],theta[1]*np.log(10),10**theta[2],theta[3]*np.log(10)]
        return np.log10(phi[0]*np.exp(-phi[1]*x)+phi[2]*np.exp(-phi[3]*x))

    def predictR(nsteps,noise,lengths,weights):
        #takes noise as sigma^2, passes on as log10
        thetagrid = nl2thetam(np.log10(noise),lengths)
        R=0
        for i in range(thetagrid.shape[1]):
            R += thetam2r(thetagrid[:,i],nsteps)*weights[i]
        return R/float(thetagrid.shape[1])

    #basic sanity check bounds
    sigmamin = para['icostfn'](Bremain) #only take 1 evaluation
    cmax = para['costfn'](sigmamin)
    acc=0
    i=0
    while acc+Oapprox(optstate.n+i)<Bremain:
        acc+=Oapprox(optstate.n+i)
        i+=1
    if i==0:
        #just go for a 1s evaluation since the overhead from this step will take us over budget
        ev['s'] = para['icostfn'](1.)
        print('end on step {} due to budget end')
        return ev['s'],persist,aux

    sigmamax = para['icostfn']((Bremain-acc)/i)
    cmin = para['costfn'](sigmamax)
    stepsmax=i
    print('maxcost {} at sigma={}, mincost {} at sigma={} for {} steps'.format(cmax,sigmamin,cmin,sigmamax,i))
    nlevels=100
    N = np.zeros([5,nlevels])
    N[0,:] = np.logspace(np.log10(sigmamin),np.log10(sigmamax),nlevels)
    N[1,:] = para['costfn'](N[0,:])
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
        N[4,i] = 10**predictR(N[2,i],N[3,i],l,A)


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
    return N[0,j],persist,aux

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
    t0=time.clock()
    ev['s'],persist,aux = getVopt(optstate,persist,para,ev,aux)
    t1=time.clock()
    print('extra due to getVopt {}'.format(t1-t0))

    return x,ev,persist,aux

for i in range(1):
    fname = 'eihyp{}.csv'.format(i)
    C=gpbo.core.config.eihypdefault(f,D,20,1e-6,fpath,fname)
    C.reccpara['kindex']=C.aqpara['kindex']=GPdc.SQUEXP
    C.reccpara['mprior']=C.aqpara['mprior']= sp.array([0.]+[0.]*D)
    C.reccpara['sprior']=C.aqpara['sprior']= sp.array([0.5]*(D+1))

    debugoutput['predictive']=False
    C.aqpara['DH_SAMPLES']=200
    C.aqpara['B']=75*cfn(1e-6)
    C.aqfn = eihyppredictiveaq
    C.stopfn = gpbo.optimize.totalTorNstopfn
    C.stoppara['nmax']=20
    C.stoppara['tmax']=C.aqpara['B']
    C.aqpara['costfn']=cfn
    C.aqpara['icostfn']=icfn


    lref,nref,M,V = pickle.load(open('/home/mark/Desktop/multidpara.p','rb'))[D]
    #mean-var of regret model
    C.aqpara['means'] = [interp2d(lref,nref,M[i].T) for i in range(len(M))]
    C.aqpara['vars']  = [interp2d(lref,nref,V[i].T) for i in range(len(V))]

    out = gpbo.search(C)
