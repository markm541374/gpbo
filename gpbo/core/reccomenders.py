# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

import scipy as sp
import os
import time
import DIRECT
import logging

logger = logging.getLogger(__name__)

import GPdc

def argminrecc(optstate,**para):
    if para['onlyafter']>len(optstate.y) or not len(optstate.y)%para['everyn']==0:
        return [sp.NaN for i in para['lb']],{'didnotrun':True}
    
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
    if para['onlyafter']>len(optstate.y) or not len(optstate.y)%para['everyn']==0:
        return [sp.NaN for i in para['lb']],{'didnotrun':True}
    logger.info('gpmap reccomender')
    d=len(para['lb'])
    x=sp.vstack(optstate.x)
    y=sp.vstack(optstate.y)
    s= sp.vstack([e['s'] for e in optstate.ev])
    dx=[e['d'] for e in optstate.ev]
    MAP = GPdc.searchMAPhyp(x,y,s,dx,para['mprior'],para['sprior'], para['kindex'])
    logger.info('MAPHYP {}'.format(MAP))
    G = GPdc.GPcore(x,y,s,dx,GPdc.kernel(para['kindex'],d,MAP))
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
                'everyn':1
                }

gpmap = gpmaprecc,gpmapprior

def gpmapasrecc(optstate,**para):
    if para['onlyafter']>len(optstate.y) or not len(optstate.y)%para['everyn']==0:
        return [sp.NaN for i in para['lb']],{'didnotrun':True}
    logger.info('gpmapas reccomender')
    d=len(para['lb'])
    
    x=sp.hstack([sp.vstack(optstate.x),sp.vstack([e['xa'] for e in optstate.ev])])
    
    
    y=sp.vstack(optstate.y)
    s= sp.vstack([e['s'] for e in optstate.ev])
    dx=[e['d'] for e in optstate.ev]
    MAP = GPdc.searchMAPhyp(x,y,s,dx,para['mprior'],para['sprior'], para['kindex'])
    logger.info('MAPHYP {}'.format(MAP))
    G = GPdc.GPcore(x,y,s,dx,GPdc.kernel(para['kindex'],d+1,MAP))
    def directwrap(xq,y):
        xq.resize([1,d])
        xe = sp.hstack([xq,sp.array([[0.]])])
        #print xe
        a = G.infer_m(xe,[[sp.NaN]])
        return (a[0,0],0)
    [xmin,ymin,ierror] = DIRECT.solve(directwrap,para['lb'],para['ub'],user_data=[], algmethod=1, volper=para['volper'], logfilename='/dev/null')
    logger.info('reccsearchresult: {}'.format([xmin,ymin,ierror]))
    
    if True:
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
                'everyn':1,
                }

gpasmap = gpmapasrecc,gpmapasprior

def gphinasrecc(optstate,**para):
    if para['onlyafter']>len(optstate.y) or not len(optstate.y)%para['everyn']==0:
        return [sp.NaN for i in para['lb']],{'didnotrun':True}
    logger.info('gpmapas reccomender')
    d=len(para['lb'])


    x=sp.hstack([sp.vstack([e['xa'] for e in optstate.ev]),sp.vstack(optstate.x)])

    
    y=sp.vstack(optstate.y)
    s= sp.vstack([e['s'] for e in optstate.ev])
    dx=[e['d'] for e in optstate.ev]
    
    G = GPdc.GPcore(x,y,s,dx,[GPdc.kernel(optstate.aux['kindex'],d+1,h) for h in optstate.aux['HYPdraws']])
    def directwrap(xq,y):
        xq.resize([1,d])
        xe = sp.hstack([sp.array([[0.]]),xq])
        #print xe
        a = G.infer_m_post(xe,[[sp.NaN]])
        return (a[0,0],0)
    [xmin,ymin,ierror] = DIRECT.solve(directwrap,para['lb'],para['ub'],user_data=[], algmethod=1, volper=para['volper'], logfilename='/dev/null')
    logger.info('reccsearchresult: {}'.format([xmin,ymin,ierror]))
    
    if True:
        if not os.path.exists(os.path.join('.', 'dbout')):
            os.mkdir('dbout')
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
            ax[0].plot(x[i,0],x[i,1],'b.')
            circle = plt.Circle([x[i,0],x[i,1]], radius=0.5*x[i,2]*l, edgecolor="none",color='lightblue',alpha=0.8-0.6*x[i,2])
            ax[0].add_patch(circle)
        ax[0].plot(x[i+1,0],x[i+1,1],'r.')
        circle = plt.Circle([x[i+1,0],x[i+1,1]], radius=0.5*x[i,2]*l, edgecolor="none",color='lightblue',alpha=0.8-0.6*x[i,2])
        ax[0].add_patch(circle)
        ax[0].axis([-1.,1.,-1.,1.])
        ax[1].plot(xmin[0],xmin[1],'ro')
        fig.savefig(os.path.join('dbout', 'datavis' + time.strftime('%d_%m_%y_%H:%M:%S') + '.png'))
        del (fig)

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
                'everyn':1,
                }

gpashin = gphinasrecc,gphinasprior