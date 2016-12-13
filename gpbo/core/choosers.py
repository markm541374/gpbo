import scipy as sp
from scipy import linalg as spl
from itertools import groupby
from collections import defaultdict
from gpbo.core import PES
import ESutils
import DIRECT
from sklearn import mixture
import logging
import gpbo
from scipy.stats import norm as norms
from gpbo.core import GPdc as GP
logger = logging.getLogger(__name__)
try:
    from matplotlib import pyplot as plt
    from matplotlib import patches
    plots=True
except ImportError:
    plots=False
    plt=None
import os
import time
from copy import deepcopy


def always0(optstate,persist,**para):
    return 0,None,dict()


def alternate(optstate,persist,**para):
    return optstate.n%2,None,dict()


def aftern(optstate,persist,**para):
    if optstate.n>para['n']:
        return 1,None,dict()
    else:
        return 0,None,dict()


def gpcommitment(optstate,persist,**para):
    logging.info('gpcommitment')
    if persist == None:
        persist = defaultdict(list)
    if optstate.n<para['onlyafter']:
        return 0,persist,dict()
    if persist['flip']:
        return 1,persist,dict()
    d = len(para['lb'])


    #build a GP with slice-samples hypers
    x = sp.vstack(optstate.x)
    y = sp.vstack(optstate.y)
    s = sp.vstack([e['s'] for e in optstate.ev])
    dx = [e['d'] for e in optstate.ev]
    G = PES.makeG(x, y, s, dx, para['kindex'], para['mprior'], para['sprior'], para['nhyp'])


    #find the pred. minimum
    def directwrap(xq,y):
        xq.resize([1,d])
        a = G.infer_m_post(xq,[[sp.NaN]])
        return (a[0,0],0)
    [xmin,ymin,ierror] = DIRECT.solve(directwrap,para['lb'],para['ub'],user_data=[], algmethod=1, maxf=para['maxf'], logfilename='/dev/null')
    optstate.startlocal = xmin

    #hess inference in rings around xmin

    #draw support points
    W = sp.vstack(ESutils.draw_support(G, para['lb'], para['ub'], para['support'], ESutils.SUPPORT_LAPAPROT, para=para['starts']))
    nd = para['draws']
    #draw from min dist over support points plus pred.min. Returns values at drawn min  and value and gradient at pred min in Y. Locations in R, indexis into W in A
    R, Y, A = ESutils.draw_min_xypairgrad(G, W, nd, xmin)


    Y_,classaux = gmmclassifier(R,xmin)
    #Y_r,classauxr = radgmmclassifier(R,xmin)
    #classify all the points in W as 0:localdrawn, 1:drawn, 2:notdrawn
    Z_ = sp.ones(para['support']+1,dtype='i')*2
    for i in xrange(nd):
        Z_[A[i]]=Y_[i]


    #regrets according to the draws
    ERegret = 0.  # full expected
    splitR = [0.,0.]  # [within,outside] the local cluster
    countnotlocal=0
    for i in xrange(nd):
        reg = max(0,Y[i, 1] - Y[i, 0])
        ERegret+=reg
        splitR[Y_[i]]+=reg
        countnotlocal+=Y_[i]

    ERegret/=float(nd)
    try:
        LRegret = splitR[0]/(float(nd))
    except ZeroDivisionError:
        LRegret=0
    try:
        GRegret = splitR[1]/(float(nd))
    except ZeroDivisionError:
        GRegret=0
    probnotlocal = countnotlocal/float(nd)

    persist['ERegret'].append(ERegret)
    persist['LRegret'].append(LRegret)
    persist['GRegret'].append(GRegret)
    persist['probnotlocal'].append(probnotlocal)

    #ER bounds on regret based on support

    ns = Z_.shape[0]
    Mf, Vf = G.infer_diag(sp.vstack([W, xmin]), [[sp.NaN]] * (W.shape[0] + 1))
    V = Vf.max(axis=0) + sp.var(Mf, axis=0).reshape([1, ns])
    M = Mf.min(axis=0)
    ER = sp.empty([1, ns])
    for i in xrange(ns):
        ER[0, i] = ERdouble(M[-1],sp.sqrt(V[0,-1]), M[i], sp.sqrt(V[0, i]))
    EBound = ER.sum()
    splits=[0,0,0]
    splitbound = [0.,0.,0.]
    splitzeros = [0,0,0]
    for i in xrange(Z_.shape[0]):
        splits[Z_[i]]+=1
        splitbound[Z_[i]]+=ER[0,i]
        if ER[0,i]==0.:
            splitzeros[Z_[i]]+=1
    GBound=splitbound[1]+splitbound[2]
    LBound=splitbound[0]

    persist['EBound'].append(EBound)
    persist['GBound'].append(GBound)
    persist['LBound'].append(LBound)

    #predictors for the bounds

    Gpred = predictforward(persist['GBound'])
    Lpred = predictforward(persist['LRegret'])
    linit = Lpred.predict(optstate.n-1)
    ginit = Gpred.predict(optstate.n-1)

    #check the switch to lacl asecision
    chpara = {
        'lrf':lambda x:min(linit,Lpred.predict(x+optstate.n-1)),
        'grf':lambda x:min(ginit,Gpred.predict(x+optstate.n-1)),
        'bcf':lambda x:0.,
        'evc':1.,
        'lsc':0.,
        'lsn':20,
        'lsr':1e-7,
        'brm':para['budget']-optstate.n,
        'tol':1e-5
        }
    now,when = choice2(chpara)
    persist['flip']=now

    from gpbo.core import debugoutput, debugoptions, debugpath
    if debugoutput and debugoptions['adaptive'] and plots:
        print 'plotting choice...'
        fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(40, 40))
        #plot the current GP
        n = 60
        x_ = sp.linspace(-1, 1, n)
        y_ = sp.linspace(-1, 1, n)
        z_ = sp.empty([n, n])
        s_ = sp.empty([n, n])
        for i in range(n):
            for j in range(n):
                m_, v_ = G.infer_diag_post(sp.array([y_[j], x_[i]]), [[sp.NaN]])
                z_[i, j] = m_[0, 0]
                s_[i, j] = sp.sqrt(v_[0, 0])

        CS = ax[0, 0].contour(x_, y_, z_, 20)
        ax[0, 0].clabel(CS, inline=1, fontsize=10)
        CS = ax[1, 0].contour(x_, y_, s_, 20)
        ax[0, 0].axis([-1., 1., -1., 1.])
        ax[1, 0].clabel(CS, inline=1, fontsize=10)
        ax[0, 0].plot(xmin[0], xmin[1], 'ro')


        #plot support points and draws
        for i in xrange(para['support']):
            symbol = ['r.','g.','bx'][Z_[i]]
            ax[0,1].plot(W[i,0],W[i,1],symbol)
            #if ER[0,i]>0.:
            #    ax[2,2].plot(W[i,0],W[i,1],'b.')
            #else:
            #    ax[2,2].plot(W[i,0],W[i,1],'r.')
        ax[0, 1].axis([-1, 1, -1, 1])

        for i in xrange(nd):
            symbol = ['r.', 'g.'][Y_[i]]
            ax[1, 1].plot(R[i, 0], R[i, 1], symbol)
        #for i in xrange(nd):
        #    symbol = ['r.', 'g.'][Y_r[i]]
        #    ax[1, 2].plot(R[i, 0], R[i, 1], symbol)

        Nselected = sorted([len(list(group)) for key, group in groupby(sorted(A))])
        Nselected.reverse()
        ax[0,2].plot(Nselected,'b')
        ax[0,2].text(1,max(Nselected)*0.5,'no draws from {} bins ({}%)'.format(para['support']-len(Nselected),100*(para['support']-len(Nselected))/float(para['support'])))

        #plot gaussian mixtures
        for i, (mean, cov) in enumerate(zip(classaux['means'], classaux['covs'])):
            color='green' if i not in classaux['localset'] else 'red'
            v, w = spl.eigh(cov)
            # Plot an ellipse to show the Gaussian component
            angle = sp.arctan2(w[0][1], w[0][0])
            angle = 180. * angle / sp.pi  # convert to degrees
            v = 2. * sp.sqrt(2.) * sp.sqrt(v)
            ell = patches.Ellipse(mean, v[0], v[1], 180. + angle, color=None, edgecolor=color,facecolor='none')
            ax[1, 1].add_artist(ell)

        #plot draw based regrets
        ax[2,0].semilogy(persist['ERegret'],'k') #full regret
        ax[2, 0].semilogy(persist['LRegret'], 'b') #local regret
        ax[2, 0].semilogy(persist['GRegret'], 'm') #nonlocal regret


        #plot support based regret bounds
        ax[2, 0].semilogy(persist['EBound'], 'k--.')  # full regret
        ax[2, 0].semilogy(persist['LBound'], 'b--.')  # local regret
        ax[2, 0].semilogy(persist['GBound'], 'm--.')  # nonlocal regret


        # display zeroER count for regret bounds
        msg = 'out of total   {}\n' \
              'total zeros    {}\n' \
              'inlocal        {}  ({})\n' \
              'drawnotlocal   {}  ({})\n' \
              'nodraw         {}  ({})\n\n' \
              'golocalnow     {}\n' \
              'switch in      {}' \
              ''.format(ns,sum(splitzeros),splitzeros[0],splits[0],splitzeros[1],splits[1],splitzeros[2],splits[2],now,when)
        ax[2, 0].text(0.5, min(persist['GBound']),msg)

        #show forward predictions
        nq=300
        xaxis = sp.linspace(0,para['budget'],nq)
        Lp = map(Lpred.predict, xaxis)
        Gp = map(Gpred.predict, xaxis)
        ax[2, 0].semilogy(xaxis,Lp, 'b:')  # local regret
        ax[2, 0].semilogy(xaxis,Gp, 'm:')


        #true regret if available:
        #try:
        Rtrue = para['cheatf'](xmin,**{'s':0.,'d':[sp.NaN]})[0]-para['cheatymin']
        persist['RTrue'].append(Rtrue)
        ax[2,1].semilogy(persist['RTrue'],'g')
        ax[2,1].semilogy(persist['ERegret'],'k') #prob not in local

        #except:
        #    pass
        try:
            fig.savefig(os.path.join(debugpath, 'lotsofplots' + time.strftime('%d_%m_%y_%H:%M:%S') + '.png'))
        except BaseException as e:
            logger.error(str(e))
        fig.clf()
        plt.close(fig)
        del (fig)

        obj = [ER,M,V,Z_,Y_,R,Y,xmin,ymin,persist,deepcopy([x, y, s, dx, para['kindex'], para['mprior'], para['sprior'], para['nhyp']])]
        import pickle
        pickle.dump(obj, open('dbout/{}.p'.format(optstate.n), 'wb'))
        logger.info('endopt')

    return int(now),persist,{'localstart':xmin}

def gmmclassifier(R,xmin):

    # fit gmm to the drawn mins
    lowest_bic = sp.inf
    bic = []
    n_components_range = range(1, 9)
    for n_components in n_components_range:
        # Fit a Gaussian mixture with EM
        gmm = mixture.GaussianMixture(n_components=n_components, covariance_type='full')
        gmm.fit(R)
        bic.append(gmm.bic(R))
        if bic[-1] < lowest_bic:
            lowest_bic = bic[-1]
            best_gmm = gmm

    clf = best_gmm

    # get the closest ellipse
    D_ = sp.inf
    local = -1
    for i, m in enumerate(clf.means_):
        dist = spl.norm(m - xmin)
        if dist < D_:
            local = i
            D_ = dist

    localset = [local]
    closestcov = clf.covariances_[local]
    closestmean = clf.means_[local]

    # add any that are similar size and overlapping
    for i, (mean, cov) in enumerate(zip(clf.means_, clf.covariances_)):
        if i not in localset:
            if (mean - closestmean).max() < 0.8*sp.sqrt(cov.max()) and (mean - closestmean).max() < 0.8*sp.sqrt(
                    closestcov.max()):
                if closestcov.max() < 1.3 * cov.max() and closestcov.max() > 0.8 * cov.max():
                    localset.append(i)

    # predict classes
    Y_ = clf.predict(R)
    #map classes to in out
    Z_ = sp.empty(Y_.size,dtype='i')
    for i in xrange(Y_.size):
        if Y_[i] in localset:
            Z_[i]=0
        else:
            Z_[i]=1
    return Z_,{'localset':localset,'covs':clf.covariances_,'means':clf.means_}

def radgmmclassifier(Ri,xmin):
    R=sp.empty([Ri.shape[0],1])
    for i in xrange(Ri.shape[0]):
        z=Ri[i,:]-xmin
        R[i,0] = sp.log10(sp.sqrt(z[0]**2+z[1]**2)) if z[0]!=0 and z[1]!=0 else -3.
    # fit gmm to the drawn mins
    lowest_bic = sp.inf
    bic = []
    n_components_range = range(1, 3)
    for n_components in n_components_range:
        # Fit a Gaussian mixture with EM
        gmm = mixture.GaussianMixture(n_components=n_components, covariance_type='full')
        gmm.fit(R)
        bic.append(gmm.bic(R))
        if bic[-1] < lowest_bic:
            lowest_bic = bic[-1]
            best_gmm = gmm

    clf = best_gmm
    local = int(sp.argmin(clf.means_))
    localset = [local]
    closestcov = clf.covariances_[local]
    closestmean = clf.means_[local]
    # add any that are similar size and overlappind

    # predict classes
    Y_ = clf.predict(R)
    #map classes to in out
    Z_ = sp.empty(Y_.size,dtype='i')
    for i in xrange(Y_.size):
        if Y_[i] in localset:
            Z_[i]=0
        else:
            Z_[i]=1

    #print "XXXXXXXXXXXXXXXX_\n{}\n{}".format(clf.means_.flatten(),clf.covariances_.flatten())
    return Z_,{'localset':localset,'covs':clf.covariances_,'means':clf.means_}

def ERdouble(m0,s0,m1,s1):
    sz = sp.sqrt(s0**2+s1**2)
    mz = m0-m1
    return sz*norms.pdf(-mz/sz)-mz*norms.cdf(-mz/sz)+mz


class predictforward:
    def __init__(self, data):
        self.n = len(data)
        self.Y = sp.array([sp.log10(i) for i in data]).reshape([self.n, 1])
        self.X = sp.array([float(i) for i in range(self.n)]).reshape([self.n, 1])
        self.D = [[sp.NaN]] * self.n
        self.S = sp.array([[1e-2]]*self.n).T
#        self.S = sp.zeros([1e-2, self.n])
#        self.ki = GP.SQUEXPCS
#        self.MAPHYP = GP.searchMAPhyp(self.X, self.Y, self.S, self.D, sp.array([1., 2., 2.]), sp.array([2., 2., 2.]),self.ki)
        self.ki = GP.CPDEC1
        self.MAPHYP = GP.searchMAPhyp(self.X, self.Y, self.S, self.D, sp.array([1., 1., 0.]), sp.array([1., 1., 0.5]),self.ki)

        self.kf = GP.kernel(self.ki, 1, sp.array(self.MAPHYP))

        self.G = GP.GPcore(self.X, self.Y, self.S, self.D, self.kf)

    def predict(self, x):
        m = self.G.infer_m(sp.array([[float(x)]]), [[sp.NaN]])
        return 10**(m[0, 0])

def choice(para):
    lrf = para['lrf'] #local regret function, in steps ahead
    grf = para['grf'] #global regret function, in steps ahead
    bcf = para['bcf'] #BayesOpt cost function, overhead cost to run BO from i to i+1 ahead
    evc = para['evc'] #Evaluation cost
    lsc = para['lsc'] #Local seach cost, overhead cost to run a local search
    lsn = para['lsn'] #Local search number, number of evaluations needed for a local search
    lsr = para['lsr'] #Local search regrert, regret (tolerance) of a local search
    brm = para['brm'] #Budget ReMaining, cost allowed to be inured

    M=0 #number of bo states available (including current)
    acc=0.
    while acc<brm:
        acc+=bcf(M)+evc
        M+=1

    Sc = sp.zeros([2,M])
    Sv = sp.zeros([2,M])
    Ac = sp.zeros(M)
    #forward pass to costs of state
    for i in xrange(1,M):
        Sc[0,i]=Sc[0,i-1]+bcf(i-1)+evc
        Sc[1,i]=Sc[0,i-1]+lsc+lsn*evc

    #backward pass for value and action
    Sv[0,M-1]=lrf(M-1)+grf(M-1)
    Sv[1,M-1]=grf(M-1)+lsr if Sc[1,M-1]<brm else sp.Inf
    switchat=-1
    for i in reversed(range(M-1)):
        Sv[1,i]=grf(i-1)+lsr if Sc[1,i]<brm else sp.Inf
        if Sv[1,i+1]>=Sv[0,i+1]:
            Sv[0,i]=Sv[0,i+1]
            Ac[i]=0
        else:
            Sv[0,i]=Sv[1,i+1]
            Ac[i]=1
            switchat=i
    #print sp.vstack([Sv,Ac]).T
    return Ac[0],switchat


def choice2(para):
    lrf = para['lrf'] #local regret function, in steps ahead
    grf = para['grf'] #global regret function, in steps ahead
    bcf = para['bcf'] #BayesOpt cost function, overhead cost to run BO from i to i+1 ahead
    evc = para['evc'] #Evaluation cost
    lsc = para['lsc'] #Local seach cost, overhead cost to run a local search
    lsn = para['lsn'] #Local search number, number of evaluations needed for a local search
    lsr = para['lsr'] #Local search regrert, regret (tolerance) of a local search
    brm = para['brm'] #Budget ReMaining, cost allowed to be inured
    tol = para['tol'] #regret tolerance
    M=0 #number of bo states available (including current)
    acc=0.
    while acc<brm:
        acc+=bcf(M)+evc
        M+=1

    Sc = sp.zeros([2,M])
    Sv = sp.zeros([2,M])
    Ac = sp.zeros(M)
    #forward pass to costs of state
    for i in xrange(1,M):
        Sc[0,i]=Sc[0,i-1]+bcf(i-1)+evc
        Sc[1,i]=Sc[0,i-1]+lsc+lsn*evc

    #backward pass for value and action
    Sv[0,M-1]=lrf(M-1)+grf(M-1)+Sc[0,M-1]
    Sv[1,M-1]=grf(M-1)+lsr +Sc[1,M-1] if Sc[1,M-1]<brm and grf(M-1)+lsr<tol else sp.Inf
    switchat=-1
    for i in reversed(range(M-1)):
        print grf(i-1)+lsr,grf(i-1)+lsr<tol,tol,Sc[1,i],Sc[1,i]<brm
        Sv[1,i]=grf(i-1)+lsr+Sc[1,i] if Sc[1,i]<brm and grf(i-1)+lsr<tol else sp.Inf
        if Sv[1,i+1]>=Sv[0,i+1]:
            Sv[0,i]=Sv[0,i+1]
            Ac[i]=0
        else:
            Sv[0,i]=Sv[1,i+1]
            Ac[i]=1
            switchat=i
    return Ac[0],switchat