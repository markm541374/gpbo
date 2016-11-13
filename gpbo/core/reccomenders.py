# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

import scipy as sp
from scipy import linalg as spl
import os
import time
import DIRECT
import itertools
import logging

logger = logging.getLogger(__name__)
try:
    from matplotlib import pyplot as plt
    from matplotlib import patches
    plots=True
except ImportError:
    plots=False
    plt=None

import GPdc

def argminrecc(optstate,**para):
    #if para['onlyafter']>len(optstate.y) or not len(optstate.y)%para['everyn']==0:
    #    return [sp.NaN for i in para['lb']],{'didnotrun':True}
    
    logger.info('argmin reccomender')
    xinc = optstate.x[0]
    yinc = 1e99
    for x,y in zip(optstate.x,optstate.y):
        if y<yinc:
            xinc = x
            yinc=y
    return xinc,{'yinc':yinc}

argminpara = dict()
argmin = argminrecc, argminpara

def gpmaprecc(optstate,**para):
    if para['onlyafter']>=len(optstate.y) or not len(optstate.y)%para['everyn']==0:
        return argminrecc(optstate,**para)
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
    [xmin,ymin,ierror] = DIRECT.solve(directwrap,para['lb'],para['ub'],user_data=[], algmethod=1, volper=para['volper'], logfilename='/dev/null')
    return [i for i in xmin],{'MAPHYP':MAP,'ymin':ymin}

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

def gpmapasrecc(optstate,**para):
    if para['onlyafter']>=len(optstate.y) or not len(optstate.y)%para['everyn']==0:
        return argminrecc(optstate, **para)
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
    [xmin,ymin,ierror] = DIRECT.solve(directwrap,para['lb'],para['ub'],user_data=[], algmethod=1, volper=para['volper'], logfilename='/dev/null')
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
        del(fig)
    return [i for i in xmin],{'MAPHYP':MAP,'ymin':ymin}

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

def gphinasrecc(optstate,**para):
    if para['onlyafter']>=len(optstate.y) or not len(optstate.y)%para['everyn']==0:
        #return [sp.NaN for i in para['lb']],{'didnotrun':True}
        return argminrecc(optstate, **para)
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
    [xmin,ymin,ierror] = DIRECT.solve(directwrap,para['lb'],para['ub'],user_data=[], algmethod=1, volper=para['volper'], logfilename='/dev/null')
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
        del(fig)


    return [i for i in xmin],{'ymin':ymin}

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


def gphinrecc(optstate,**para):
    print [para['onlyafter'],len(optstate.y)]
    if para['onlyafter']>=len(optstate.y) or not len(optstate.y)%para['everyn']==0:
        #return [sp.NaN for i in para['lb']],{'didnotrun':True}
        return argminrecc(optstate, **para)
    logger.info('gphin reccomender')
    d=len(para['lb'])


    x=sp.vstack(optstate.x)

    
    y=sp.vstack(optstate.y)
    s= sp.vstack([e['s'] for e in optstate.ev])
    dx=[e['d'] for e in optstate.ev]
    print optstate.aux
    G = GPdc.GPcore(x, y, s, dx, [GPdc.kernel(optstate.aux['kindex'], d, h) for h in optstate.aux['HYPdraws']])
    def directwrap(xq,y):
        xq.resize([1,d])
        xe = xq
        #print xe
        a = G.infer_m_post(xe,[[sp.NaN]])
        return (a[0,0],0)
    [xmin,ymin,ierror] = DIRECT.solve(directwrap,para['lb'],para['ub'],user_data=[], algmethod=1, volper=para['volper'], logfilename='/dev/null')
    logger.info('reccsearchresult: {}'.format([xmin,ymin,ierror]))

    
    return [i for i in xmin],{'ymin':ymin}

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



def adaptiverecc(optstate,**para):
    if para['onlyafter']>=len(optstate.y) or not len(optstate.y)%para['everyn']==0:
        return argminrecc(optstate,**para)
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
    [xmin,ymin,ierror] = DIRECT.solve(directwrap,para['lb'],para['ub'],user_data=[], algmethod=1, volper=para['volper'], logfilename='/dev/null')


    import ESutils
    W = sp.vstack(ESutils.draw_support(G, para['lb'], para['ub'], para['support'], ESutils.SUPPORT_LAPAPROT, para=para['starts']))
    nd = para['draws']

    R,Y = ESutils.draw_min_xypairgrad(G, W, nd, xmin)

    A = sp.empty([nd,1])
    for i in xrange(nd):
        A[i,0] = -Y[i,2:].dot(R[i,:]-xmin)/(spl.norm(Y[i,2:])*spl.norm(R[i,:]-xmin))
        #print '__________________\n'
        #print xmin
        #print Y[i,2:]
        #print Y[i,2:].dot(xmin)
        #print [spl.norm(Y[i,2:]),spl.norm(xmin)]
        #print A[i,0]

    D = sp.zeros([nd,1])
    for i in xrange(nd):
        for j in xrange(len(para['ub'])):
            D[i,0]+=(xmin[j]-R[i,j])**2
        D[i,0]=sp.sqrt(D[i,0])
    from gpbo.core import debugoutput
    from gpbo.core import debugoptions
    if debugoutput and debugoptions['support'] and plots:
        print 'plotting support...'
        fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(30, 20))

        n = 100
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

        ax[0,0].plot(xmin[0],xmin[1],'ro')

        ax[0,1].plot(W[:, 0], W[:, 1], 'bx')
        ax[0,1].plot(R[:, 0], R[:, 1], 'r.')


        axA = ax[1, 1].twinx()
        lowD = 1e2
        for i in xrange(nd):
            if 0.<D[i,0]<lowD and A[i,0]<0:
                lowD=D[i,0]
            axA.plot(D[i,0],A[i,0],'g.')
        axA.plot([lowD,lowD],[-1,1],'g')

        zerocount = 0
        ERegret = 0.
        IRegret = 0.
        ORegret = 0.
        r = [0, 0]

        for i in xrange(nd):
            if D[i,0]==0.:
                zerocount+=1
            ax[1,1].plot(D[i,0],Y[i,0],'b.')
            ax[1, 1].plot(D[i, 0], Y[i, 1], 'r.')

            ERegret+=Y[i,1]-Y[i,0]
            if D[i,0]<lowD:
                IRegret+=Y[i,1]-Y[i,0]
                r[0]+=1
            else:
                ORegret+=Y[i,1]-Y[i,0]
                r[1]+=1

        ERegret/=float(nd)
        try:
            IRegret/=float(r[0])
        except ZeroDivisionError:
            IRegret=-1
        try:
            ORegret /= float(r[1])
        except ZeroDivisionError:
            Oregret=-1
        prob = r[0]/float(nd)
        ax[1,1].plot([0,D.max()],[ymin,ymin],'r')



        ax[1,1].set_xscale('log')
        ax[1, 1].text(ax[1,1].get_xlim()[0], ymin, str(zerocount))
        #ax[1, 1].text(ax[1, 1].get_xlim()[0]*2., ymin, str(ERegret))

        try:
            CRegret = para['cheatf'](xmin,**{'s': 0., 'd': [sp.NaN]})[0]-para['cheatymin']
        except KeyError:
            CRegret=sp.nan

        msg = 'ExpR full   : {}\n' \
              'ExpR in loc : {}\n' \
              'ExpR out loc: {}\n' \
              'Pmin in loc : {}\n' \
              'TrueR       : {}'.format(ERegret,IRegret,ORegret,prob,CRegret)
        l_,u_=ax[1, 1].get_xlim()
        axA.text(l_+0.02*(u_-l_),-0.98,msg)

        #---------------------------------------------------
        """
        from sklearn.cluster import MeanShift, estimate_bandwidth
        bandwidth = estimate_bandwidth(R, quantile=0.2, n_samples=0.25*R.shape[0])
        ms = MeanShift(bandwidth=bandwidth, bin_seeding=True, cluster_all=False)
        ms.fit(R)
        labels = ms.labels_
        cluster_centers = ms.cluster_centers_

        labels_unique = sp.unique(labels)
        n_clusters_ = len(labels_unique)
        if labels_unique[0] == -1:
            n_clusters_ -= 1

        print("number of estimated clusters : %d" % n_clusters_)

        from itertools import cycle
        colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
        for k, col in zip(range(-1, n_clusters_), colors):
            my_members = labels == k
            if k >= 0:
                ax[0,2].plot(R[my_members, 0], R[my_members, 1], col + '.')
            else:
                ax[0,2].plot(R[my_members, 0], R[my_members, 1], '.', color='grey')
        """
        #---------------------------------------------------
        from sklearn import mixture
        lowest_bic = sp.inf
        bic = []
        n_components_range = range(1, 7)
        cv_types = ['spherical', 'tied', 'diag', 'full']
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

        D = sp.inf
        local = -1
        for i,m in enumerate(clf.means_):
            dist=spl.norm(m-xmin)
            if dist<D:
                local=i
                D=dist

        splot=ax[0,2]

        Y_ = clf.predict(R)
        for i, (mean, cov) in enumerate(zip(clf.means_, clf.covariances_)):
            color='blue' if i != local else 'red'
            if typ == 'spherical':
                cov = sp.diag([cov] * R.shape[1])
            elif typ == 'tied' or typ == 'diag':
                cov = sp.diag(cov)
            else:
                pass
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

        from gpbo.core import debugpath
        fig.savefig(os.path.join(debugpath, 'support' + time.strftime('%d_%m_%y_%H:%M:%S') + '.png'))

        del (fig)
    return [i for i in xmin],{'MAPHYP':MAP,'ymin':ymin}
