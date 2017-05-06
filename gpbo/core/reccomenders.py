# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.
from __future__ import print_function
xrange=range
import scipy as sp
import os
import sys
import copy
from scipy.optimize import minimize
from gpbo.core.optutils import multilocal
import time
import gpbo
from gpbo.core.optutils import silentdirect as direct
import logging
logger = logging.getLogger(__name__)
try:
    from matplotlib import pyplot as plt
    from matplotlib import patches
    plots=True
except ImportError:
    plots=False
    plt=None

from gpbo.core import GPdc

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
    lb = para['lb']
    ub = para['ub']
    maxf = para['maxf']
    x=sp.vstack(optstate.x)
    y=sp.vstack(optstate.y)
    s= sp.vstack([e['s']+10**optstate.condition for e in optstate.ev])
    dx=[e['d'] for e in optstate.ev]
    MAP = GPdc.searchMAPhyp(x, y, s, dx, para['mprior'], para['sprior'], para['kindex'])
    logger.info('MAPHYP {}'.format(MAP))
    G = GPdc.GPcore(x, y, s, dx, GPdc.kernel(para['kindex'], d, MAP))
#
#    if para['smode'] == 'direct':
#        def directwrap(xq, y):
#            xq.resize([1, d])
#            a = G.infer_m_post(xq,[[sp.NaN]])
#            return (a[0, 0], 0)
#
#        [xmin, ymin, ierror] = direct(directwrap, lb, ub, user_data=[], algmethod=1, maxf=para['maxf'], logfilename='/dev/null')
#        logger.info('DIRECT found min post at {} {} {}'.format(xmin, ymin, ierror))
#    elif para['smode'] == 'multi':
#        def multiwrap(x):
#            xq = copy.copy(x)
#            xq.resize([1, d])
#            a = G.infer_m_post(xq,[[sp.NaN]])
#            return a[0, 0]
#
#        [xmin, ymin, ierror] = multilocal(multiwrap, lb, ub, maxf=maxf)
#
#        logger.info('multilocal found min post at {} {} {}'.format(xmin, ymin, ierror))
#    elif para['smode'] == 'dthenl':
#        def directwrap(xq, y):
#            xq.resize([1, d])
#            a = G.infer_m_post(xq,[[sp.NaN]])
#            return (a[0, 0], 0)
#
#        [dxmin, dymin, ierror] = direct(directwrap, lb, ub, user_data=[], algmethod=1, maxf=maxf - 150,
#                                        logfilename='/dev/null')
#        logger.info('DIRECT found min post at {} {} {}'.format(dxmin, dymin, ierror))
#
#        def localwrap(x):
#            xq = copy.copy(x)
#            xq.resize([1, d])
#            a = G.infer_m_post(xq,[[sp.NaN]])
#            return a[0, 0]
#
#        res = minimize(localwrap, dxmin, method='L-BFGS-B', bounds=tuple([(lb[j], ub[j]) for j in range(d)]),
#                       options={'ftol': 0.00001, 'maxfun': 150})
#        xmin, ymin, ierror = res.x, res.fun, res.message
#        logger.info('localrefine found min post at {} {} {}'.format(xmin, ymin, ierror))
#
#    else:
#        raise KeyError('not a search mode')
    def wrap(x):
        xq = copy.copy(x)
        xq.resize([1, d])
        a = G.infer_m_post(xq,[[sp.NaN]])
        return a[0, 0]
    xmin,ymin,ierror = gpbo.core.optutils.twopartopt(wrap,para['lb'],para['ub'],para['dpara'],para['lpara'])
    return [i for i in xmin],persist,{'MAPHYP':MAP,'ymin':ymin}

def gpmap2upperrecc(optstate,persist,**para):
    if para['onlyafter']>=len(optstate.y) :
        print('{} <= {} : switch to argmin'.format(len(optstate.y),para['onlyafter']))
        return argminrecc(optstate,persist,**para)
        #return [sp.NaN for i in para['lb']],{'didnotrun':True}
    logger.info('gpmapucb2 reccomender')
    d=len(para['lb'])
    x=sp.vstack(optstate.x)
    y=sp.vstack(optstate.y)
    s= sp.vstack([e['s']+10**optstate.condition for e in optstate.ev])
    dx=[e['d'] for e in optstate.ev]
    MAP = GPdc.searchMAPhyp(x, y, s, dx, para['mprior'], para['sprior'], para['kindex'])
    logger.info('MAPHYP {}'.format(MAP))
    G = GPdc.GPcore(x, y, s, dx, GPdc.kernel(para['kindex'], d, MAP))
    count=0
    def wrap(xq):
        xq.resize([1,d])
        a,v = G.infer_diag_post(xq,[[sp.NaN]])
        return a[0,0]+2.*sp.sqrt(v[0,0])

    #print('nevals={}\n\n'.format(count))
    #[xmin,ymin,ierror] = direct(directwrap,para['lb'],para['ub'],user_data=[], algmethod=1, maxf=para['maxf'], logfilename='/dev/null')
    xmin,ymin,ierror = gpbo.core.optutils.twopartopt(wrap,para['lb'],para['ub'],para['dpara'],para['lpara'])
    logger.info('DIRECT found post. min {} at {} {}'.format(ymin,xmin,ierror))
    a,v= G.infer_diag_post(sp.array([[xmin]]),[[sp.NaN]])
    print('mean {} std{} '.format(a[0,0],sp.sqrt(v[0,0])))

    x,p,di = argminrecc(optstate,persist,**para)
    a, v = G.infer_diag_post(sp.array([[x]]), [[sp.NaN]])
    print('at argmin mean {} std {} wasactually {}'.format(a[0, 0], sp.sqrt(v[0, 0]),di['yinc']))
    import sys
    sys.stdout.flush()
    return [i for i in xmin],persist,{'MAPHYP':MAP,'ymin':ymin}

def gpmapasrecc(optstate,persist,**para):
    if para['onlyafter']>=len(optstate.y) or not len(optstate.y)%para['everyn']==0:
        return argminrecc(optstate,persist, **para)
        #return [sp.NaN for i in para['lb']],{'didnotrun':True}
    logger.info('gpmapas reccomender')
    d=len(para['lb'])
    
    x=sp.hstack([sp.vstack(optstate.x),sp.vstack([e['xa'] for e in optstate.ev])])
    
    
    y=sp.vstack(optstate.y)
    s= sp.vstack([e['s']+10**optstate.condition for e in optstate.ev])
    dx=[e['d'] for e in optstate.ev]
    MAP = GPdc.searchMAPhyp(x, y, s, dx, para['mprior'], para['sprior'], para['kindex'])
    logger.info('MAPHYP {}'.format(MAP))
    G = GPdc.GPcore(x, y, s, dx, GPdc.kernel(para['kindex'], d + 1, MAP))
    def wrap(xq):
        xq.resize([1,d])
        xe = sp.hstack([xq,sp.array([[0.]])])
        a = G.infer_m_post(xe,[[sp.NaN]])
        return a[0,0]
    #[xmin,ymin,ierror] = direct(directwrap,para['lb'],para['ub'],user_data=[], algmethod=1, maxf=para['maxf'], logfilename='/dev/null')
    xmin,ymin,ierror = gpbo.core.optutils.twopartopt(wrap,para['lb'],para['ub'],para['dpara'],para['lpara'])
    logger.info('reccsearchresult: {}'.format([xmin,ymin,ierror]))


    from gpbo.core import debugoutput
    if debugoutput['datavis']:
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

def gphinasargminrecc(optstate, persist, **para):
    if para['onlyafter'] >= len(optstate.y) or not len(optstate.y) % para['everyn'] == 0:
        # return [sp.NaN for i in para['lb']],{'didnotrun':True}
        return argminrecc(optstate, persist, **para)
    logger.info('gpmapas reccomender')
    d = len(para['lb'])

    x = sp.hstack([sp.vstack([e['xa'] for e in optstate.ev]), sp.vstack(optstate.x)])

    y = sp.vstack(optstate.y)
    s = sp.vstack([e['s']+10**optstate.condition for e in optstate.ev])
    dx = [e['d'] for e in optstate.ev]

    G = GPdc.GPcore(x, y, s, dx, [GPdc.kernel(optstate.aux['kindex'], d + 1, h) for h in optstate.aux['HYPdraws']])

    def wrap(xq):
        xq.resize([1, d])
        xe = sp.hstack([sp.array([[0.]]), xq])
        # print xe
        a = G.infer_m_post(xe, [[sp.NaN]])
        return a[0, 0]
    best=sp.Inf
    incumbent = None
    for i in range(len(optstate.x)):
        thisone = wrap(sp.array(optstate.x[i]))
        if thisone<best:
            best=thisone
            incumbent=optstate.x[i]

    logger.info('reccsearchresult: x {} pred.y {}'.format(incumbent,best))



    return [i for i in incumbent], persist, {'ymin': best}


def gphinasrecc(optstate,persist,**para):
    if para['onlyafter']>=len(optstate.y) or not len(optstate.y)%para['everyn']==0:
        #return [sp.NaN for i in para['lb']],{'didnotrun':True}
        return argminrecc(optstate,persist, **para)
    logger.info('gpmapas reccomender')
    d=len(para['lb'])

    x=sp.hstack([sp.vstack([e['xa'] for e in optstate.ev]),sp.vstack(optstate.x)])
    y=sp.vstack(optstate.y)
    s= sp.vstack([e['s']+10**optstate.condition for e in optstate.ev])
    dx=[e['d'] for e in optstate.ev]
    
    G = GPdc.GPcore(x, y, s, dx, [GPdc.kernel(optstate.aux['kindex'], d + 1, h) for h in optstate.aux['HYPdraws']])
#    def directwrap(xq,y):
#        xq.resize([1,d])
#        xe = sp.hstack([sp.array([[0.]]),xq])
#        #print xe
#        a = G.infer_m_post(xe,[[sp.NaN]])
#        return (a[0,0],0)
#    [xmin,ymin,ierror] = direct(directwrap,para['lb'],para['ub'],user_data=[], algmethod=1, maxf=para['maxf'], logfilename='/dev/null')

    def wrap(x):
        xq = copy.copy(x)
        xq.resize([1, d])
        xe = sp.hstack([sp.array([[0.]]),xq])
        a = G.infer_m_post(xe,[[sp.NaN]])
        return a[0, 0]

    xmin,ymin,ierror = gpbo.core.optutils.twopartopt(wrap,para['lb'],para['ub'],para['dpara'],para['lpara'])
    logger.info('reccsearchresult: {}'.format([xmin,ymin,ierror]))
    from gpbo.core import debugoutput
    if debugoutput['datavis']:
        if not os.path.exists(debugoutput['path']):
            os.mkdir(debugoutput['path'])

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
        fig.savefig(os.path.join(debugoutput['path'], 'datavis' + time.strftime('%d_%m_%y_%H:%M:%S') + '.png'))
        fig.clf()
        plt.close(fig)
        del(fig)


    return [i for i in xmin],persist,{'ymin':ymin}

def gphinrecc(optstate,persist,**para):
    #print( [para['onlyafter'],len(optstate.y)])
    if para['onlyafter']>=len(optstate.y) or not len(optstate.y)%para['everyn']==0:
        #return [sp.NaN for i in para['lb']],{'didnotrun':True}
        return argminrecc(optstate, persist,**para)
    logger.info('gphin reccomender')
    d=len(para['lb'])
    x=sp.vstack(optstate.x)

    y=sp.vstack(optstate.y)
    s= sp.vstack([e['s']+10**optstate.condition for e in optstate.ev])
    dx=[e['d'] for e in optstate.ev]
   # print optstate.aux
    G = GPdc.GPcore(x, y, s, dx, [GPdc.kernel(optstate.aux['kindex'], d, h) for h in optstate.aux['HYPdraws']])
    lb = para['lb']
    ub = para['ub']

    def wrap(x):
        xq = copy.copy(x)
        xq.resize([1, d])
        a = G.infer_m_post(xq,[[sp.NaN]])
        return a[0, 0]

    xmin,ymin,ierror = gpbo.core.optutils.twopartopt(wrap,lb,ub,para['dpara'],para['lpara'])
    return [i for i in xmin],persist,{'ymin':ymin}

