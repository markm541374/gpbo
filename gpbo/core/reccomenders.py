# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.
from __future__ import print_function
import scipy as sp
from scipy import linalg as spl
from scipy import stats as sps
import os
import sys
import time

from OPTutils import silentdirect as direct
import itertools
import logging
from scipy.optimize import minimize as spm
logger = logging.getLogger(__name__)
try:
    from matplotlib import pyplot as plt
    from matplotlib import patches
    plots=True
except ImportError:
    plots=False
    plt=None

import GPdc

def argminrecc(optstate,persist,**para):
    #if para['onlyafter']>len(optstate.y) or not len(optstate.y)%para['everyn']==0:
    #    return [sp.NaN for i in para['lb']],{'didnotrun':True}
    
    logger.info('argmin reccomender')
    xinc = optstate.x[0]
    yinc =sp.Inf
    for x,y in zip(optstate.x,optstate.y):
        if y<yinc:
            xinc = x
            yinc=y
    return xinc,persist,{'yinc':yinc}

argminpara = dict()
argmin = argminrecc, argminpara

def gpmaprecc(optstate,persist,**para):
    if para['onlyafter']>=len(optstate.y) or not len(optstate.y)%para['everyn']==0:
        return argminrecc(optstate,persist,**para)
        #return [sp.NaN for i in para['lb']],{'didnotrun':True}
    logger.info('gpmap reccomender')
    d=len(para['lb'])
    x=sp.vstack(optstate.x)
    y=sp.vstack(optstate.y)
    s= sp.vstack([e['s'] for e in optstate.ev])
    dx=[e['d'] for e in optstate.ev]
    MAP = GPdc.searchMAPhyp(x, y, s, dx, para['mprior'], para['sprior'], para['kindex'])
    logger.info('MAPHYP {}'.format(MAP))
    G = GPdc.GPcore(x, y, s, dx, GPdc.kernel(para['kindex'], d, MAP))
    def directwrap(xq,y):
        xq.resize([1,d])
        a = G.infer_m(xq,[[sp.NaN]])
        return (a[0,0],0)


    [xmin,ymin,ierror] = direct(directwrap,para['lb'],para['ub'],user_data=[], algmethod=1, volper=para['volper'], logfilename='/dev/null')


    logger.info('DIRECT found post. min {} at {} {}'.format(ymin,xmin,ierror))
    return [i for i in xmin],persist,{'MAPHYP':MAP,'ymin':ymin}

gpmapprior = {
                'ev':{'s':1e-9,'d':[sp.NaN]},
                'lb':[-1.-1.],
                'ub':[1.,1.],
                'mprior':sp.array([1.,0.,0.]),
                'sprior':sp.array([1.,1.,1.]),
                'kindex':GPdc.MAT52,
                'volper':1e-6,
                'onlyafter':10,
                'check':False,
                'everyn':1
                }

gpmap = gpmaprecc,gpmapprior

def gpmapasrecc(optstate,persist,**para):
    if para['onlyafter']>=len(optstate.y) or not len(optstate.y)%para['everyn']==0:
        return argminrecc(optstate,persist, **para)
        #return [sp.NaN for i in para['lb']],{'didnotrun':True}
    logger.info('gpmapas reccomender')
    d=len(para['lb'])
    
    x=sp.hstack([sp.vstack(optstate.x),sp.vstack([e['xa'] for e in optstate.ev])])
    
    
    y=sp.vstack(optstate.y)
    s= sp.vstack([e['s'] for e in optstate.ev])
    dx=[e['d'] for e in optstate.ev]
    MAP = GPdc.searchMAPhyp(x, y, s, dx, para['mprior'], para['sprior'], para['kindex'])
    logger.info('MAPHYP {}'.format(MAP))
    G = GPdc.GPcore(x, y, s, dx, GPdc.kernel(para['kindex'], d + 1, MAP))
    def directwrap(xq,y):
        xq.resize([1,d])
        xe = sp.hstack([xq,sp.array([[0.]])])
        #print xe
        a = G.infer_m(xe,[[sp.NaN]])
        return (a[0,0],0)
    [xmin,ymin,ierror] = direct(directwrap,para['lb'],para['ub'],user_data=[], algmethod=1, volper=para['volper'], logfilename='/dev/null')
    logger.info('reccsearchresult: {}'.format([xmin,ymin,ierror]))
    from gpbo.core import debugoutput, debugoptions
    if debugoutput and debugoptions['datavis']:
        A2 = MAP[0]
        l = MAP[3]
        from matplotlib import pyplot as plt
        fig, ax = plt.subplots( nrows=3, ncols=1 ,figsize=(10,30))
        
        
        n = 200
        x_ = sp.linspace(-1,1,n)
        y_ = sp.linspace(-1,1,n)
        z_ = sp.empty([n,n])
        s_ = sp.empty([n,n])
        for i in xrange(n):
            for j in xrange(n):
                m_,v_ = G.infer_diag_post(sp.array([y_[j],x_[i],0.]),[[sp.NaN]])
                z_[i,j] = m_[0,0]
                s_[i,j] = sp.sqrt(v_[0,0])
        CS = ax[1].contour(x_,y_,z_,20)
        ax[1].clabel(CS, inline=1, fontsize=10)
        CS = ax[2].contour(x_,y_,s_,20)
        ax[2].clabel(CS, inline=1, fontsize=10)
        for i in xrange(x.shape[0]):
            ax[0].plot(x[i,0],x[i,1],'b.')
            circle = plt.Circle([x[i,0],x[i,1]], radius=0.5*x[i,2]*l, edgecolor="none",color='lightblue',alpha=0.8-0.6*x[i,2])
            ax[0].add_patch(circle)
            
        ax[0].axis([-1.,1.,-1.,1.])
        fig.savefig(os.path.join(os.path.expanduser('~'),'Dropbox/debugoutput','datavis'+time.strftime('%d_%m_%y_%H:%M:%S')+'.png'))
        fig.clf()
        plt.close(fig)
        del(fig)
    return [i for i in xmin],persist,{'MAPHYP':MAP,'ymin':ymin}

gpmapasprior = {
                'ev':{'s':1e-9,'d':[sp.NaN],'xa':0.},
                #'ev':[1e-9,[sp.NaN]],
                'lb':[-1.-1.],
                'ub':[1.,1.],
                'mprior':sp.array([1.,0.,0.,0.]),
                'sprior':sp.array([1.,1.,1.,1.]),
                'kindex':GPdc.MAT52,
                'volper':1e-6,
                'onlyafter':10,
                'check':False,
                'everyn':1,
                }

gpasmap = gpmapasrecc,gpmapasprior

def gphinasrecc(optstate,persist,**para):
    if para['onlyafter']>=len(optstate.y) or not len(optstate.y)%para['everyn']==0:
        #return [sp.NaN for i in para['lb']],{'didnotrun':True}
        return argminrecc(optstate,persist, **para)
    logger.info('gpmapas reccomender')
    d=len(para['lb'])


    x=sp.hstack([sp.vstack([e['xa'] for e in optstate.ev]),sp.vstack(optstate.x)])

    
    y=sp.vstack(optstate.y)
    s= sp.vstack([e['s'] for e in optstate.ev])
    dx=[e['d'] for e in optstate.ev]
    
    G = GPdc.GPcore(x, y, s, dx, [GPdc.kernel(optstate.aux['kindex'], d + 1, h) for h in optstate.aux['HYPdraws']])
    def directwrap(xq,y):
        xq.resize([1,d])
        xe = sp.hstack([sp.array([[0.]]),xq])
        #print xe
        a = G.infer_m_post(xe,[[sp.NaN]])
        return (a[0,0],0)
    [xmin,ymin,ierror] = direct(directwrap,para['lb'],para['ub'],user_data=[], algmethod=1, volper=para['volper'], logfilename='/dev/null')
    logger.info('reccsearchresult: {}'.format([xmin,ymin,ierror]))

    from gpbo.core import debugoutput, debugoptions
    if debugoutput and debugoptions['datavis']:
        from gpbo.core import debugpath
        if not os.path.exists(debugpath):
            os.mkdir(debugpath)

        l = sp.mean([h[3] for h in optstate.aux['HYPdraws']])
        from matplotlib import pyplot as plt
        fig, ax = plt.subplots( nrows=3, ncols=1 ,figsize=(10,30))
        
        
        n = 200
        x_ = sp.linspace(-1,1,n)
        y_ = sp.linspace(-1,1,n)
        z_ = sp.empty([n,n])
        s_ = sp.empty([n,n])
        for i in xrange(n):
            for j in xrange(n):
                m_,v_ = G.infer_diag_post(sp.array([0.,y_[j],x_[i]]),[[sp.NaN]])
                z_[i,j] = m_[0,0]
                s_[i,j] = sp.sqrt(v_[0,0])
        CS = ax[1].contour(x_,y_,z_,20)
        ax[1].clabel(CS, inline=1, fontsize=10)
        CS = ax[2].contour(x_,y_,s_,20)
        ax[2].clabel(CS, inline=1, fontsize=10)
        for i in xrange(x.shape[0]-1):
            ax[0].plot(x[i,1],x[i,2],'b.')
            circle = plt.Circle([x[i,1],x[i,2]], radius=0.5*x[i,0]*l, edgecolor="none",color='lightblue',alpha=0.8-0.6*x[i,2])
            ax[0].add_patch(circle)
        ax[0].plot(x[i+1,1],x[i+1,2],'r.')
        circle = plt.Circle([x[i+1,1],x[i+1,2]], radius=0.5*x[i,0]*l, edgecolor="none",color='lightblue',alpha=0.8-0.6*x[i,2])
        ax[0].add_patch(circle)
        ax[0].axis([-1.,1.,-1.,1.])
        ax[1].plot(xmin[0],xmin[1],'ro')
        fig.savefig(os.path.join(debugpath, 'datavis' + time.strftime('%d_%m_%y_%H:%M:%S') + '.png'))
        fig.clf()
        plt.close(fig)
        del(fig)


    return [i for i in xmin],persist,{'ymin':ymin}

gphinasprior = {
                'ev':{'s':1e-9,'d':[sp.NaN],'xa':0.},
                #'ev':[1e-9,[sp.NaN]],
                'lb':[-1.-1.],
                'ub':[1.,1.],
                #'mprior':sp.array([1.,0.,0.,0.]),
                #'sprior':sp.array([1.,1.,1.,1.]),
                #'kindex':GPdc.MAT52,
                'volper':1e-6,
                'onlyafter':10,
                'check':False,
                'everyn':1,
                }

gpashin = gphinasrecc,gphinasprior


def gphinrecc(optstate,persist,**para):
    print( [para['onlyafter'],len(optstate.y)])
    if para['onlyafter']>=len(optstate.y) or not len(optstate.y)%para['everyn']==0:
        #return [sp.NaN for i in para['lb']],{'didnotrun':True}
        return argminrecc(optstate, persist,**para)
    logger.info('gphin reccomender')
    d=len(para['lb'])


    x=sp.vstack(optstate.x)

    
    y=sp.vstack(optstate.y)
    s= sp.vstack([e['s'] for e in optstate.ev])
    dx=[e['d'] for e in optstate.ev]
   # print optstate.aux
    G = GPdc.GPcore(x, y, s, dx, [GPdc.kernel(optstate.aux['kindex'], d, h) for h in optstate.aux['HYPdraws']])
    def directwrap(xq,y):
        xq.resize([1,d])
        xe = xq
        #print xe
        a = G.infer_m_post(xe,[[sp.NaN]])
        return (a[0,0],0)


    [xmin,ymin,ierror] = direct(directwrap,para['lb'],para['ub'],user_data=[], algmethod=1, volper=para['volper'], logfilename='/dev/null')


    logger.info('DIRECT found post. min {} at {} {}'.format(ymin,xmin,ierror))

    return [i for i in xmin],persist,{'ymin':ymin}

gphinprior = {
                'ev':{'s':1e-9,'d':[sp.NaN],'xa':0.},
                #'ev':[1e-9,[sp.NaN]],
                'lb':[-1.-1.],
                'ub':[1.,1.],
                #'mprior':sp.array([1.,0.,0.,0.]),
                #'sprior':sp.array([1.,1.,1.,1.]),
                #'kindex':GPdc.MAT52,
                'volper':1e-5,
                'onlyafter':10,
                'check':False,
                'everyn':1,
                }

gphin = gphinrecc,gphinprior



def adaptiverecc(optstate,persist,**para):
    if persist==None:
        persist={'ERCom':[],'probcluster':[],'Rtrue':[],'RCtrue':[],'ER':[],'ERlocal':[],'globRbound':[],'globCbound':[],'LT':[],'globRupper':[],'globCupper':[],'sumER':[]}
    if para['onlyafter']>=len(optstate.y) or not len(optstate.y)%para['everyn']==0:
        return argminrecc(optstate,persist,**para)
        #return [sp.NaN for i in para['lb']],{'didnotrun':True}
    logger.info('gpmap reccomender')
    d=len(para['lb'])
    x=sp.vstack(optstate.x)
    y=sp.vstack(optstate.y)
    s= sp.vstack([e['s'] for e in optstate.ev])
    dx=[e['d'] for e in optstate.ev]
    MAP = GPdc.searchMAPhyp(x, y, s, dx, para['mprior'], para['sprior'], para['kindex'])
    logger.info('MAPHYP {}'.format(MAP))
    #G = GPdc.GPcore(x, y, s, dx, GPdc.kernel(para['kindex'], d, MAP))
    from gpbo.core import PES
    G=PES.makeG(x, y, s, dx,para['kindex'],sp.array([1.]+[0.]*d),sp.array([1.]*(d+1)),8)

    def directwrap(xq,y):
        xq.resize([1,d])
        a = G.infer_m(xq,[[sp.NaN]])
        return (a[0,0],0)
    [xmin,ymin,ierror] = direct(directwrap,para['lb'],para['ub'],user_data=[], algmethod=1, volper=para['volper'], logfilename='/dev/null')

    def spowrap(x):
        a = G.infer_m(x,[[sp.NaN]])[0,0]
        #z = obj(x,0.,[sp.NaN])
        return a
    res = spm(spowrap, xmin,  method='l-bfgs-b',bounds=[(-1,1),(-1,1)],options={'ftol':1e-10})
    xmin = res['x']
    ymin = res['fun']

    def directwrap(xq,y):
        xq.resize([1,d])
        e = G.infer_EI_post(xq,[[sp.NaN]],wrt=ymin)
        return (-e[0,0],0)
    [xcmax,cmax,ierror] = direct(directwrap,para['lb'],para['ub'],user_data=[], algmethod=1, volper=para['volper'], logfilename='/dev/null')

    def spowrap(x):
        e = G.infer_EI_post(x,[[sp.NaN]])
        #z = obj(x,0.,[sp.NaN])
        return -e[0,0]
    res = spm(spowrap, xcmax,  method='l-bfgs-b',bounds=[(-1,1),(-1,1)],options={'ftol':1e-10})
    xcmax = res['x']
    cmax = -res['fun']
    yatcmax = G.infer_m(xcmax,[[sp.NaN]])[0,0]

    import ESutils
    W = sp.vstack(ESutils.draw_support(G, para['lb'], para['ub'], para['support'], ESutils.SUPPORT_LAPAPROT, para=para['starts']))
    nd = para['draws']

    R,Y,A = ESutils.draw_min_xypairgrad(G, W, nd, xmin)

    M,V=G.infer_diag_post(W,[[sp.NaN]]*W.shape[0])
    ER = G.infer_EI_post(W,[[sp.NaN]]*W.shape[0],wrt=ymin)
    ERc = G.infer_EI_post(W,[[sp.NaN]]*W.shape[0],wrt=ymin)
    #A = sp.empty([nd,1])
    #for i in xrange(nd):
    #    A[i,0] = -Y[i,2:].dot(R[i,:]-xmin)/(spl.norm(Y[i,2:])*spl.norm(R[i,:]-xmin))

    D = sp.zeros([nd,1])
    for i in xrange(nd):
        for j in xrange(len(para['ub'])):
            D[i,0]+=(xmin[j]-R[i,j])**2
        D[i,0]=sp.sqrt(D[i,0])
    from gpbo.core import debugoutput
    from gpbo.core import debugoptions
    if debugoutput and debugoptions['adaptive'] and plots:
        print( 'plotting support...')
        fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(40, 40))

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

        CS = ax[0,0].contour(x_, y_, z_, 20)
        ax[0,0].clabel(CS, inline=1, fontsize=10)
        CS = ax[1,0].contour(x_, y_, s_, 20)
        ax[0,0].axis([-1., 1., -1., 1.])
        ax[1,0].clabel(CS, inline=1, fontsize=10)


        ax[1,0].plot(xcmax[0],xcmax[1],'bo')
        ax[0,0].plot(xmin[0],xmin[1],'ro')

        ax[0,1].plot(W[:, 0], W[:, 1], 'bx')
        ax[0,1].plot(R[:, 0], R[:, 1], 'r.')

        #---------------------------------------------------
        from sklearn import mixture
        lowest_bic = sp.inf
        bic = []
        n_components_range = range(1, 7)
        cv_types = ['full']
        for cv_type in cv_types:
            for n_components in n_components_range:
                # Fit a Gaussian mixture with EM
                gmm = mixture.GaussianMixture(n_components=n_components,
                                              covariance_type=cv_type)
                gmm.fit(R)
                bic.append(gmm.bic(R))
                if bic[-1] < lowest_bic:
                    lowest_bic = bic[-1]
                    typ = cv_type
                    best_gmm = gmm

        clf = best_gmm

        D_ = sp.inf
        local = -1
        for i,m in enumerate(clf.means_):
            dist=spl.norm(m-xmin)
            if dist<D_:
                local=i
                D_=dist

        localset = [local]
        closestcov = clf.covariances_[local]

        closestmean = clf.means_[local]
        for i, (mean, cov) in enumerate(zip(clf.means_, clf.covariances_)):
            if i not in localset:
                if (mean - closestmean).max() < sp.sqrt(cov.max()) and (mean - closestmean).max()<sp.sqrt(closestcov.max()):
                    if closestcov.max() < 1.5 * cov.max() and closestcov.max() > 0.5 * cov.max():
                            localset.append(i)

        splot=ax[0,2]

        Y_ = clf.predict(R)
        W_ = clf.predict(W)
        for i, (mean, cov) in enumerate(zip(clf.means_, clf.covariances_)):
            color='blue' if i not in localset else 'red'
            v, w = spl.eigh(cov)
            if not sp.any(Y_ == i):
                continue
            splot.scatter(R[Y_ == i, 0], R[Y_ == i, 1], .8, color=color)

            # Plot an ellipse to show the Gaussian component
            angle = sp.arctan2(w[0][1], w[0][0])
            angle = 180. * angle / sp.pi  # convert to degrees
            v = 2. * sp.sqrt(2.) * sp.sqrt(v)
            ell = patches.Ellipse(mean, v[0], v[1], 180. + angle, color=None, edgecolor=color,facecolor='none')
            splot.add_artist(ell)

        ERegret = 0. #full expected
        IRegret = 0. #within the local cluster
        ORegret = 0. #outside the cluster



        StoConv=[]
        #NRegret = 0. #within the cluster and not converging
        r = [0, 0, 0]
        for i in xrange(nd):
            if Y_[i] in localset:

                ax[1, 2].plot(D[i, 0], Y[i, 1] - Y[i, 0], 'r.')
                IRegret+=Y[i, 1]-Y[i, 0]
                r[0]+=1

                # zero the in Y_ elements of ER
                for j in xrange(para['support']):
                    if all([R[i,k]==W[j,k] for k in range(d)]):
                        ER[0,j]=0.

            else:

                ax[1, 2].plot(D[i, 0], Y[i, 1] - Y[i, 0], 'b.')
                ORegret += Y[i, 1] - Y[i, 0]
                r[1]+=1
            ERegret += Y[i, 1] - Y[i, 0]

        ERegret /= float(nd)
        try:
            IRegret /= float(r[0])
        except ZeroDivisionError:
            IRegret = 0
        try:
            ORegret /= float(r[1])
        except ZeroDivisionError:
            ORegret = 0

        #try:
        #    NRegret /= float(r[2])
        #except ZeroDivisionError:
        #    NRegret = 0
        Iprob = r[0] / float(nd)
        #Nprob = r[2] / float(r[0])
        try:
            CRegret = para['cheatf'](xmin, **{'s': 0., 'd': [sp.NaN]})[0] - para['cheatymin']
            LocalTrace=[CRegret]
            def f(x):
                y=para['cheatf'](x, **{'s': 0., 'd': [sp.NaN]})[0]
                LocalTrace.append(y-para['cheatymin'])
                return y
            from scipy.optimize import minimize
            resob = minimize(f,xmin,method='l-bfgs-b',bounds=[(-1,1),(-1,1)],options={'ftol':1e-10})
            RCtrue = (resob['fun']-para['cheatymin'],resob['nfev'])

            print( 'optstart {} {}\n' \
                  'optend   {} {}'.format(xmin,ymin,resob['x'],resob['fun']))
        except KeyError:
            CRegret = sp.nan
            RCtrue = (sp.nan,0,)
        CommitRegret=(1-Iprob)*ORegret


        Rbound = ER.sum()

        #panyless = 1.- (1 - c_) ** float(para['support'])
        #panyless_2terms = para['support'] * c_ - 0.5 * para['support'] * (para['support'] - 1) * c_ ** 2
        #if panyless == 0.:
        #    Rupper = panyless_2terms * (Rbound / (c_))
        #else:
        #    Rupper = panyless * (Rbound / (c_))

        approxcommitregret = max(CommitRegret, Rbound)


        #uppermsg = 'XXXXXXXX\n' \
        #           'cdf upto ymin      : {}\n' \
        #           'nsupp points       : {}\n' \
        #           'ER                 : {}\n' \
        #           'ER | lessthan      : {}\n' \
        #           'p allgreaterthan   : {}\n' \
        #           'p2termless         : {}\n' \
        #           'p anylessthan      : {}\n' \
        #           'ER after           : {}\n' \
        #           ''.format(c_,para['support'],Rbound,Rbound/(c_),(1-c_)**float(para['support']),panyless_2terms,panyless,Rupper)

        #print "{} {} {} {} {}".format(S_,c_,p_,cmax,yatcmax)
        persist['RCtrue'].append(RCtrue)
        persist['ER'].append(ERegret)
        persist['ERCom'].append(CommitRegret)
        persist['probcluster'].append(Iprob)
        persist['Rtrue'].append(CRegret)
        persist['ERlocal'].append(IRegret)
        persist['globRbound'].append(Rbound)
        #persist['globRupper'].append(Rupper)
        #persist['globCbound'].append(1.-c_)
        #persist['globCupper'].append(panyless if panyless>0. else panyless_2terms)
        persist['LT'].append(LocalTrace)
        persist['sumER'].append(ER.sum())
        ax[2,0].semilogy(persist['ER'],'k')
        ax[2, 1].semilogy(persist['ER'], 'k')
        ax[2, 1].semilogy(persist['Rtrue'],'g')
        for i in xrange(len(persist['Rtrue'])):
            #ax[2,1].plot([i,i+persist['RCtrue'][i][1]],[persist['Rtrue'][i],persist['RCtrue'][i][0]],'g')
            ax[2, 1].plot([i,len(persist['LT'][i])+i], [persist['LT'][i][0],persist['LT'][i][-1]], 'r')
        #ax[2, 0].semilogy(persist['globRupper'], 'orange')

        ax[2,0].semilogy(persist['ERlocal'],'b')
        #print 'XXXXXXXXX{}'.format(persist['globRbound'])
        ax[2,0].semilogy(persist['globRbound'],'pink')
        ax[2, 0].semilogy(persist['ERCom'], 'm.-')

        #ax[2, 0].semilogy(persist['sumER'], 'cyan')
        ax[2,2].semilogy([1-i for i in persist['globCbound']],'m')
        ax[2,2].semilogy([1-i for i in persist['probcluster']])
        ax[2, 2].semilogy(persist['globCupper'],'orange')

        ax[2,2].set_ylabel('probcluster')

        l=len(persist['Rtrue'])
        #mguess = 0.01 / (xmin.size + 1)
        #mx=0
        #for i in xrange(nd):
        #    if Y_[i] == local:
        #        stepsahead=(-sp.log(Y[i, 1] - Y[i, 0]) - 10 * sp.log(10)) / sp.log(mguess)
        #        if mx<stepsahead:
        #            mx=stepsahead


        ax[2,0].plot([l-1,l+20-1],[ERegret,approxcommitregret],'grey')
        ax[2, 0].plot([l - 1, l + 20 - 1], [approxcommitregret, approxcommitregret], 'grey')
                #ax[2, 0].plot([l + stepsahead - 1], [approxcommitregret], 'ko')
        #ax[2, 0].plot([l-1, l-1 + persist['RCtrue'][l-1][1]], [persist['Rtrue'][l-1], persist['RCtrue'][l-1][0]], 'g')
        ax[2, 0].plot(range(l-1, len(LocalTrace)+l-1), LocalTrace, 'g')
        msg = 'TrueR        : {}\n' \
              'ExpR full    : {}\n\n' \
              'Prob in C    : {}\n' \
              'ExpR in C    : {}\n' \
              'ExpR notin C : {}\n\n' \
              'ComR         : {}\n' \
              ''.format(CRegret, ERegret, Iprob, IRegret, ORegret, CommitRegret)
        l_x, u_x = ax[1, 1].get_xlim()
        l_y, u_y = ax[1, 1].get_ylim()
        ax[1,1].text(l_x + 0.02 * (u_x - l_x), l_y + 0.02 * (u_y - l_y), '\n----\n'+msg,fontdict={'size':18})
        # ---------------------------------------------------

        from gpbo.core import debugpath
        fig.savefig(os.path.join(debugpath, 'lotsofplots' + time.strftime('%d_%m_%y_%H:%M:%S') + '.png'))
        fig.clf()
        plt.close(fig)
        del (fig)

        import pickle
        obj = [W,R,Y,ERc,Y_,localset,CommitRegret,IRegret,ORegret,ERegret,optstate.n,persist,clf.means_, clf.covariances_,M,V,xmin,ymin,W_,A]
        pickle.dump(obj, open('dbout/{}.p'.format(optstate.n), 'wb'))
        logger.info('endopt')
        print( 'endstep')


    return [i for i in xmin],persist,{'MAPHYP':MAP,'ymin':ymin}
