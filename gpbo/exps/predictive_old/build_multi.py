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
parser.add_argument('-m', '--mode', dest='mode', action='store', default=2,type=int)

args = parser.parse_args()

D = 2#args.dimension
lb = np.array([-1.]*D)
ub = np.array([1.]*D)
#lengthscales from 0.05 to 1.5
#lengthscale = [np.round(gammadist.rvs(3.,scale=0.15),3) for i in range(D)]
lengthscale = [0.2]*D
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
#fpath = ['all','twice','once','none'][args.mode]
from gpbo.core import debugoutput


gammalogpdf = lambda x,shape,scale: (shape-1)*np.log(x)-x/scale#-shape*np.log(scale)-np.log(gamma(shape))
normlogpdf = lambda x,loc,scale: -((x-loc)/scale)**2
lognormlogpdf = lambda x,s,loc,scale: -((np.log(x)-loc)/scale)**2 - np.log(x)
tlogpdf = lambda x,v,scale: -0.5*(v+1)*np.log(1+((x/scale)**2)/v)

def ccmodel(x,theta):
    y0 = theta[1]
    g0 = theta[2]
    y1 = theta[3]
    g1 = theta[4]
    a=0.
    b=200.
    r=200.
    r2 = 40000.
    r3 = 8000000.
    c0 = 3*y0/r2 + g0/r
    d0 = 2*y0/r3 + g0/r2
    c1 = 3*y1/r2 - g1/r
    d1 = -2*y1/r3 + g1/r2
    Y = ((theta[5]-1.)*np.exp(-theta[6]*x)+1.)*((d0*(x-b)+c0)*(x-b)**2 + (d1*(x-a)+c1)*(x-a)**2)
    return np.select([x<10,x>=10],[theta[0],Y])

def extraT(X,Y,P):
    theta0 = np.empty(10)
    for i in range(10):
        if i in [1,2,3,4,5,6,7,8,9]:
            theta0[i] = P[i][0]*P[i][2]
        elif i in [0]:
            theta0[i] = -P[i][1]
        else:
            theta0[i] = P[i][0]

    if Y.size == 0:
        return theta0
    if Y.size < 10:
        thetaopt = theta0
        thetaopt[0] = np.mean(Y)
        thetaopt[9] = np.sqrt(np.var(Y))
        return thetaopt

    #print(theta0)
    def llk(theta):
        pa=0
        for i in range(10):
            if i in [0,1,2,3,4,5,6,7,8,9]:
                #print(gammalogpdf(theta[i],P[i][0],P[i][2]))
                pa += gammalogpdf(theta[i],P[i][0],P[i][2])
            else:
                pa += lognormlogpdf(theta[i],P[i][0],P[i][1],P[i][2])
        Yh = ccmodel(X,theta)
        E = Y-Yh
        L=-pa
        L-=np.sum(normlogpdf(E[:10],0.,theta[9]))
        L-=np.sum(tlogpdf(E[10:],1./theta[7],theta[8]))
        #L-=np.sum(normlogpdf(E[10:],0.,theta[8]))
        return L
    #print(llk(theta0))
    bds = ((1e-9,np.Inf),(1e-9,np.Inf),(1e-9,np.Inf),(1e-9,np.Inf),
              (1e-9,np.Inf),(1e-9,np.Inf),(1e-9,np.Inf),(1e-9,np.Inf),
             (1e-9,np.Inf),(1e-9,np.Inf))
    res = minimize(llk,theta0,method='L-BFGS-B',bounds=bds)
    #print(res)
    return res.x

def overlearn(X,Y,P):
    Theta = extraT(X,Y,P)
    f = lambda x:ccmodel(x,Theta)
    return f


def getVopt2(optstate,persist,para,ev,aux, init=False):
    P = pickle.load(open('overheadv3.p','rb'))[0]
    D = len(para['lb'])
    #overhead predictor
    opredict = overlearn(np.arange(optstate.n),np.array(optstate.aqtime),P)
    OIapprox = lambda x:np.sum([opredict(i) for i in range(x)])
    Oapprox = lambda x:opredict(x)
    #hyperparamteres
    if init==False:
        A = np.array([h[0] for h in aux['HYPdraws']])
        L = np.empty([D,para['DH_SAMPLES']])
        for i in range(D):
            L[i,:] = [h[i+1] for h in aux['HYPdraws']]
        l_ = np.power(np.prod(L,axis=0),1./float(D))
        if D!=2:
            raise #next line is 2d
        l_ = np.array([np.sqrt(h[1]*h[2]) for h in aux['HYPdraws']])

    else:
        A = np.ones(600)
        if para['priorshape']=='gamma':
            L = np.empty([D,600])
            for i in range(D):
                L[i,:] = stats.gamma.rvs(para['mprior'][i+1], size=600, scale=para['sprior'][i+1], loc=0.)
            l_ = np.power(np.prod(L,axis=0),1./float(D))

        else:
            raise
    ind = np.argsort(l_)
    l=l_[ind]
    A=A[ind]
    ldist = gammadist(*gammadist.fit(l_,floc=0))

    means = para['means']
    vars = para['vars']
    nhyp = A.size#len(thetam[0])
    Bgone = sum(optstate.aqtime)+optstate.C
    B = para['B']
    Bremain = B-Bgone

    def nl2thetam(n,l):
        return np.array([p(l,n) for p in means])

    def thetam2r(theta,x):
        phi = [10**theta[0],theta[1]*np.log(10),10**theta[2],theta[3]*np.log(10)]
        return np.log10(phi[0]*np.exp(-phi[1]*x)+phi[2]*np.exp(-phi[3]*x))

    def predictR(nsteps,noise,ldist):
        lengths = np.linspace(0.,ldist.ppf(1-1e-6),400)
        h=lengths[1]-lengths[0]
        weights = ldist.pdf(lengths)
        #takes noise as sigma^2, passes on as log10
        thetagrid = nl2thetam(np.log10(noise),lengths)
        R=0
        for i in range(thetagrid.shape[1]):
            R += 10**(thetam2r(thetagrid[:,i],nsteps))*weights[i]*h
        return R
    steps = np.arange(1,200-optstate.n)
    over = np.cumsum([Oapprox(s+optstate.n) for s in steps])
    BperS = (Bremain-over)/steps
    nlevels = np.argmax(BperS<0)
    sigmaval = icfn(BperS)

    if nlevels==0:
        #just go for a 1s evaluation since the overhead from this step will take us over budget
        ev['s'] = para['icostfn'](1.)
        print('end on step {} due to budget end')
        return ev['s'],persist,aux

    print('maxcost {} at sigma={}, mincost {} at sigma={} for {} steps'.format(BperS[0],sigmaval[0],BperS[nlevels-1],sigmaval[nlevels-1],nlevels+optstate.n))
    N = np.zeros([5,nlevels])
    N[0,:] = sigmaval[:nlevels]
    N[1,:] = BperS[:nlevels]
    N[2,:]=steps[:nlevels]+optstate.n

    s = np.sum([np.log(e['s']) for e in optstate.ev])
    #use the logmean variance under projected policy as noise level to predict performance
    N[3,:] = np.exp(((N[2,:]-optstate.n)*np.log(N[0,:]) + s)/N[2,:])
    for i in range(nlevels):
        N[4,i] = predictR(N[2,i],N[3,i],ldist)


    j = np.argmin(N[4,:])
    ev['s'] = N[0,j]
    expectedsteps = N[2,j]
    expectedcost = N[1,j]
    if init:
        lastaq=np.NaN
    else:
        lastaq=optstate.aqtime[-1]
    print('expected final regret {} after {} steps. Next ov prediction {} evcost {} meanl {} meanA {} lastover {}'.format(N[4,j],expectedsteps,Oapprox(optstate.n),expectedcost,np.mean(l),np.mean(A),lastaq))

    aux['lmean']=np.mean(l)
    aux['steppredict']=expectedsteps
    aux['overpredict']=Oapprox(optstate.n)
    aux['regretpredict']=N[4,j]

    if debugoutput['predictive']:
        from matplotlib import pyplot as plt
        f,a = plt.subplots(nrows=5,ncols=2,figsize=[14,8])
        #a[0,0].plot(np.sort(A),np.linspace(0,1,A.size))
        if not 'stepexpect' in persist.keys():
            persist['stepexpect']=np.array([[optstate.n,N[2,j]]])
        else:
            persist['stepexpect'] = np.vstack([persist['stepexpect'],[optstate.n,N[2,j]]])
        z = np.linspace(0.,ldist.ppf(1-1e-6),400)
        a[4,0].hist(L[0,:],20,normed=1)
        dist = gammadist(*gammadist.fit(L[0,:],floc=0))
        a[4,0].plot(z,dist.pdf(z),'r')
        a[4,1].hist(L[1,:],20,normed=1)
        dist = gammadist(*gammadist.fit(L[1,:],floc=0))
        a[4,1].plot(z,dist.pdf(z),'r')

        a[0,0].hist(l,20,normed=1)

        dist = gammadist(*gammadist.fit(l,floc=0))
        a[0,0].plot(z,dist.pdf(z),'r')
        a[0,1].hist(A,20,normed=1)
        a[0,1].set_ylabel('A')
        a[0,0].set_ylabel('l')
        a[1,1].semilogx(N[0,:],N[2,:],'b.-')
        a[1,1].plot(N[0,j],N[2,j],'ro')
        a[1,0].loglog(N[0,:],N[4,:],'b.')
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
    persist['s']=N[0,j]
    return N[0,j],persist,aux
def eihyppredictiveaq(optstate,persist,**para):
    x,ev,persist,aux = gpbo.core.acquisitions.eihypaq(optstate,persist,**para)
    t0=time.clock()
    if para['vchecks']=='none':
        ev['s']=para['v']
        return x,ev,persist,aux
    elif optstate.n==0:
        #S = np.empty([40,2])
        #for i in range(40):
        #   S[i,0] =  getVopt2(optstate,persist,para,ev,aux,init=True)[0]
        #   S[i,1] =  getVopt(optstate,persist,para,ev,aux,init=True)[0]
        #print(S)
        #print(np.mean(S,axis=0),np.sqrt(np.var(S,axis=0)))
        #raise KeyError
        ev['s'],persist,aux = getVopt2(optstate,persist,para,ev,aux,init=True)
        t1=time.clock()
        print('extra due to getVopt {}'.format(t1-t0))
        return x,ev,persist,aux
    elif (optstate.n==para['nrandinit'] and para['vchecks']=='twice') or (para['vchecks']=='all' and optstate.n>=para['nrandinit']):
        ev['s'],persist,aux = getVopt2(optstate,persist,para,ev,aux)
        t1=time.clock()
        print('extra due to getVopt {}'.format(t1-t0))
        return x,ev,persist,aux
    else:
        persist['overhead']=time.clock()-t0
        ev['s'] = persist['s']
        return x,ev,persist,aux
    return

for i in range(1):
    fname = 'eihyp{}_l_{}_{}.csv'.format(args.index, int(1000 * lengthscale[0]), int(1000 * lengthscale[1]))
    fname = 'eihyp{}.csv'.format(i)
    C=gpbo.core.config.eihypdefault(f,D,20,1e-6,fpath,fname)
    C.reccpara['kindex']=C.aqpara['kindex']= GPdc.SQUEXP

    C.reccpara['mprior']=C.aqpara['mprior']= sp.array([2.]+[3.]*D)
    C.reccpara['sprior']=C.aqpara['sprior']= sp.array([0.5]+[0.15]*D)
    C.reccpara['priorshape']=C.aqpara['priorshape']='gamma'
    C.reccpara['onlyafter']=C.aqpara['nrandinit']= 10
    debugoutput['predictive']=True
    C.aqpara['DH_SAMPLES']=100
    C.aqpara['B']=100*cfn(1e-6)
    C.aqfn = eihyppredictiveaq
    C.stopfn = gpbo.optimize.totalTorNstopfn
    C.stoppara['nmax']=200
    C.stoppara['tmax']=C.aqpara['B']
    C.aqpara['costfn']=cfn
    C.aqpara['icostfn']=icfn

    C.aqpara['vchecks']=['all','twice','once','none'][0]#[args.mode]
    if C.aqpara['vchecks']=='none':
        ntarget = np.random.randint(50,150)
        C.aqpara['v']=icfn(C.aqpara['B']/float(ntarget))
    lref,nref,M,V = pickle.load(open('multidpara.p','rb'))[D]
    #mean-var of regret model
    C.aqpara['means'] = [interp2d(lref,nref,M[i].T) for i in range(len(M))]
    C.aqpara['vars']  = [interp2d(lref,nref,V[i].T) for i in range(len(V))]

    out = gpbo.search(C)
