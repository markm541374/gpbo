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


def always0(optstate,persist,**para):
    return 0


def alternate(optstate,persist,**para):
    return optstate.n%2


def aftern(optstate,persist,**para):
    if optstate.n>para['n']:
        return 1
    else:
        return 0


def gpcommitment(optstate,persist,**para):
    logging.info('gpcommitment')
    if persist == None:
        persist = defaultdict(list)
    if optstate.n<para['onlyafter']:
        return 0,persist
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


    #draw support points
    W = sp.vstack(ESutils.draw_support(G, para['lb'], para['ub'], para['support'], ESutils.SUPPORT_LAPAPROT, para=para['starts']))
    nd = para['draws']
    #draw from min dist over support points plus pred.min. Returns values at drawn min  and value and gradient at pred min in Y. Locations in R, indexis into W in A
    R, Y, A = ESutils.draw_min_xypairgrad(G, W, nd, xmin)


    Y_,classaux = gmmclassifier(R,xmin)
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

    ER = G.infer_EI_post(sp.vstack([W,xmin]), [[sp.NaN]] * (W.shape[0]+1), wrt=ymin)
    EBound = ER.sum()

    splitbound = [0.,0.,0.]
    splitzeros = [0,0,0]
    for i in xrange(Z_.shape[0]):
        splitbound[Z_[i]]+=ER[0,i]
        if ER[0,i]==0.:
            splitzeros[Z_[i]]+=1
    GBound=splitbound[1]+splitbound[2]
    LBound=splitbound[0]

    persist['EBound'].append(EBound)
    persist['GBound'].append(GBound)
    persist['LBound'].append(LBound)


    #more conservative? bounds
    ns = Z_.shape[0]
    Mf, Vf = G.infer_diag(sp.vstack([W, xmin]), [[sp.NaN]] * (W.shape[0] + 1))
    V2 = Vf.max(axis=0) + sp.var(Mf, axis=0).reshape([1, ns])
    M2 = Mf.min(axis=0)
    ER2 = sp.empty([1, ns])
    for i in xrange(ns):
        ER2[0, i] = gpbo.core.GPdc.EI(ymin, M2[i], V2[0, i])

    EBound2 = ER2.sum()

    splitbound2 = [0., 0., 0.]
    splitzeros2 = [0, 0, 0]
    closetoV = [0,0,0]
    for i in xrange(Z_.shape[0]):
        splitbound2[Z_[i]] += ER2[0, i]
        if ER2[0, i] == 0.:
            splitzeros2[Z_[i]] += 1
   
    GBound2 = splitbound2[1] + splitbound2[2]
    LBound2 = splitbound2[0]

    persist['EBound2'].append(EBound2)
    persist['GBound2'].append(GBound2)
    persist['LBound2'].append(LBound2)

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
        ax[0, 1].axis([-1, 1, -1, 1])

        for i in xrange(nd):
            symbol = ['r.', 'g.'][Y_[i]]
            ax[1, 1].plot(R[i, 0], R[i, 1], symbol)

        Nselected = sorted([len(list(group)) for key, group in groupby(sorted(A))])
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
        ax[2,1].semilogy(persist['probnotlocal'],'k') #prob not in local


        #plot support based regret bounds
        ax[2, 0].semilogy(persist['EBound'], 'k:x')  # full regret
        ax[2, 0].semilogy(persist['LBound'], 'b:x')  # local regret
        ax[2, 0].semilogy(persist['GBound'], 'm:x')  # nonlocal regret

        #display zeroER count for regret bounds

        ax[2,0].text(0.5,min(persist['GBound']),'out of total  {}\ntotal zeros {}\ninlocal {}\ndrawnnotlocal {}\nnotdrawn {}'.format(Z_.shape[0],sum(splitzeros),splitzeros[0],splitzeros[1],splitzeros[2]))

        #repeat for more conservative
        # plot support based regret bounds
        ax[2, 0].semilogy(persist['EBound2'], 'k:o')  # full regret
        ax[2, 0].semilogy(persist['LBound2'], 'b:o')  # local regret
        ax[2, 0].semilogy(persist['GBound2'], 'm:o')  # nonlocal regret

        # display zeroER count for regret bounds

        ax[2, 0].text(0.5*len(persist['GBound2'])+0.5, min(persist['GBound']),'out of total  {}\ntotal zeros {}\ninlocal {}\ndrawnnotlocal {}\nnotdrawn {}'.format(Z_.shape[0],sum(splitzeros2),splitzeros2[0],splitzeros2[1],splitzeros2[2]))


        fig.savefig(os.path.join(debugpath, 'lotsofplots' + time.strftime('%d_%m_%y_%H:%M:%S') + '.png'))
        fig.clf()
        plt.close(fig)
        del (fig)

        ER = G.infer_EI_post(sp.vstack([W, xmin]), [[sp.NaN]] * (W.shape[0] + 1), wrt=ymin)
        M,V = G.infer_diag_post(sp.vstack([W, xmin]), [[sp.NaN]] * (W.shape[0] + 1))

        obj = [ER,M,V,Z_,Y_,R,Y,ymin,persist,Mf,Vf]
        import pickle
        pickle.dump(obj, open('dbout/{}.p'.format(optstate.n), 'wb'))
        logger.info('endopt')

    return 0,persist

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