from __future__ import print_function
import scipy as sp
from scipy import linalg as spl
from scipy.optimize import minimize
from itertools import groupby
from collections import defaultdict
from gpbo.core import PES
from gpbo.core import ESutils
import DIRECT
from gpbo.core.optutils import silentdirectwrapped as direct
from sklearn import mixture
import logging
import tqdm
import gpbo
from scipy.stats import norm as norms
from gpbo.core import GPdc as GP
logger = logging.getLogger(__name__)
try:
    from matplotlib import pyplot as plt
    from matplotlib import patches
    plots=True
    plt.style.use('seaborn-paper')
except ImportError:
    plots=False
    plt=None
import os
import time
from copy import deepcopy


def always0(optstate,persist,**para):
    return 0,None,dict()

def globallocalregret(optstate,persist,**para):
    try:
        return globallocalregret_(optstate,persist,**para)
    except gpbo.core.GPdc.MJMError as e:
        if persist['raiseS']:
            persist['raiseS']+=1
        else:
            persist['raiseS']=-20
        logger.error('error in globallocalregret, inflating noise to {}\n{}'.format(persist['raiseS'],e))
        return globallocalregret_(optstate,persist,**para)

def globallocalregret_(optstate,persist,**para):
    #doublenormdist
    #norprior
    if persist == None:
        persist = defaultdict(list)
        persist['raiseS']=False
    if optstate.n < para['onlyafter']:
        return 0, persist, dict()
    if persist['flip']:
        return 1, persist, dict()

    logging.info('globallocalregretchooser with {} inflated diagonal'.format(persist['raiseS']))
    d = len(para['lb'])
    lb = para['lb']
    ub = para['ub']

    #build a GP with slice-samples hypers
    x = sp.vstack(optstate.x)
    y = sp.vstack(optstate.y)
    if persist['raiseS']:
        s = sp.vstack([e['s']+10**persist['raiseS'] for e in optstate.ev])
    else:
        s = sp.vstack([e['s'] for e in optstate.ev])
    dx = [e['d'] for e in optstate.ev]
    logger.info('building GP')
    G = PES.makeG(x, y, s, dx, para['kindex'], para['mprior'], para['sprior'], para['nhyp'])


    #find the pred. minimum

    xmin,ymin,ierror = gpbo.core.optutils.twopartopt(lambda x:G.infer_m_post(x,[[sp.NaN]])[0,0],para['lb'],para['ub'],para['dpara'],para['lpara'])
    logger.info('post min {} at {} ({})'.format(xmin,ymin,ierror))
    xvmax,vmax,ierror = gpbo.core.optutils.twopartopt(lambda x:-G.infer_diag_post(x,[[sp.NaN]])[1][0,0],para['lb'],para['ub'],para['dpara'],para['lpara'])
    mvmax,vvmax = [j[0,0] for j in G.infer_diag_post(xvmax,[[sp.NaN]])]
    logger.info('post var max {} at {} with mean {} ({})'.format(vvmax,xvmax,mvmax,ierror))

    #get hessian/grad posterior
    #local probmin elipse at post min
    GH = gpbo.core.optutils.gpGH(G,xmin)
    Gr,cG,H,Hvec,varHvec,M,varM = GH

    #est the local regret
    Mdraws = gpbo.core.GPdc.draw(M[0,:],varM,200)
    lrest=0.
    for i in xrange(200):
        sM = Mdraws[i,:]
        sG = sM[:d]
        sH = gpbo.core.optutils.Hvec2H(sM[d:],d)
        sR = 0.5*sG.dot(sp.linalg.solve(sH,sG))
        lrest+= max(0.,sR)
    lrest/=200.
    logger.info('localregretest {}'.format(lrest))

    #step out to check +ve defininteness
    Rad = sp.logspace(para['pveballrrange'][0],para['pveballrrange'][1],para['pveballrsteps'])
    PP = -sp.ones(Rad.size)
    logger.info('checking for +ve definite ball')
    pc = gpbo.core.optutils.probgppve(G,sp.array(xmin),tol=para['pvetol'])
    logger.info('prob pvedef at xmin {}'.format(pc))

    PDcondition = lambda x:gpbo.core.optutils.probgppve(G,sp.array(x)+sp.array(xmin),tol=para['pvetol'])>1-para['pvetol']

    rmax = gpbo.core.optutils.ballradsearch(d,1.,PDcondition,neval=100,lineSmax=20)

    logger.info('+ve region radius {}'.format(rmax))
    if gpbo.core.debugoptions['adaptive']:
        fig, ax = plt.subplots(nrows=3, ncols=4, figsize=(85, 85))
        # plot the current GP
        if d==2:
            gpbo.core.optutils.gpplot(ax[0,0],ax[0,1],G,para['lb'],para['ub'],ns=60)
            ax[0,0].set_title('GP_post_mean')
            ax[0,1].set_title('GP_post_var')
            ax[0, 0].plot(xmin[0], xmin[1], 'ro')
            #plot some draws from H
            for i in xrange(20):
                Gm,Gv,Hd = gpbo.core.drawconditionH(*GH)
                try:
                    sp.linalg.cholesky(Hd)
                    gpbo.core.optutils.plotprobstatellipse(Gv,Hd,xmin,ax[1,1],logr=True)
                except sp.linalg.LinAlgError:
                    pass
        if rmax>0:
            ax[1,1].plot([sp.log10(rmax)]*2,[0.,2*sp.pi],'purple')
        else:
            logger.debug('plotting some draws...')
            #draw support points
            W = sp.vstack(ESutils.draw_support(G, lb, ub, 2000, ESutils.SUPPORT_LAPAPROT, para=20))
            nd = 1500
            #draw mins and value of g at xmin as pair
            R, Y, A = ESutils.draw_min_xypairgrad(G, W, nd, xmin)
            #plot support
            if d==2:
                gpbo.core.optutils.plotaslogrtheta(W[:,0],W[:,1],xmin[0],xmin[1],ax[1,1],'b.')
                ax[0,2].plot(W[:,0],W[:,1],'b.')
                #plot mindraws
                gpbo.core.optutils.plotaslogrtheta(R[:,0],R[:,1],xmin[0],xmin[1],ax[1,1],'r.')
                ax[0,2].plot(R[:,0],R[:,1],'r.')
        ax[1,3].semilogx(Rad,PP)
    if rmax==0:
        if gpbo.core.debugoptions['adaptive']:
            try:
                from gpbo.core import debugpath
                fname = 'lotsofplots' + time.strftime('%d_%m_%y_%H:%M:%S') + '.png'
                print('saving as {}'.format(fname))
                fig.savefig(os.path.join(debugpath, fname))
            except BaseException as e:
                logger.error(str(e))
            fig.clf()
            plt.close(fig)
            del (fig)
        logger.info('no +ve def region, choosereturns 0')
        return 0,persist,{'reuseH':[k.hyp for k in G.kf]}
    #draw support points
    W = sp.vstack(ESutils.draw_support(G, lb, ub, para['support'], ESutils.SUPPORT_LAPAPROT, para=20))

    Q = gpbo.core.optutils.drawpartitionmin(G,W,xmin,rmax,para['draws'])

    #pcurves from Q
    def data2cdf(X):
        n = X.size
        C = sp.linspace(1./n,1,n)
        XC = sorted(X)
        return XC,C

    Yin,Cin = data2cdf(Q[:,1])
    normin = sp.stats.norm.fit(Yin)

    Yat,Cat = data2cdf(Q[:,3])
    normat = sp.stats.norm.fit(Yat)

    Yout,Cout = data2cdf(Q[:,2])

    #normal dist with var same as max in gp model and passing through estimated prob of min sample
    ydrawmin=Yout[0]
    cdfymin=Cout[0]
    mu = ydrawmin-sp.sqrt(vvmax*2)*sp.special.erfinv(2*cdfymin-1.)
    logger.info('upper norm at y {} c {} has mu {},var {}'.format(ymin,cdfymin,mu,vvmax))
    logger.info('lower norm at x {} has mu {},var {}'.format(xvmax,mvmax,vvmax))

    #interpolator for cdf
    def splicecdf(y):
        if y<Yout[0]:
            return sp.stats.norm.cdf(y,loc=mu,scale=sp.sqrt(vvmax))
        elif y>=Yout[-1]:
            return 1.-1e-20
        else:
            i=0
            while Yout[i]<y:
                i+=1
            return Cout[i]
        return


    binormdist = gpbo.core.optutils.bigaussmin(sp.mean(Yout),sp.sqrt(sp.var(Yout)),mvmax,sp.sqrt(vvmax),0.)
    binormdistF = gpbo.core.optutils.bigaussmin(*gpbo.core.optutils.bigaussmin().fit(sp.array(Yout)))

    racc = 0.
    m,v=normin
    n=len(Cout)
    nstd=para['tailnstd'] # number of std dev below min sample to integrate over
    tailsupport = para['tailsupport']
    #regret from samples after the min
    for i in xrange(1,n):
        racc+= gpbo.core.GPdc.EI(-Yout[i],-m,v)[0,0]/float(n)
    tmp=racc
    #regret from the tail bound
    for i,y in enumerate(sp.linspace(Yout[0]-nstd*sp.sqrt(vvmax),Yout[0],tailsupport)):
        racc+= gpbo.core.GPdc.EI(-y,-m,v)[0,0]*sp.stats.norm.pdf(y,mu,sp.sqrt(vvmax))*(nstd*sp.sqrt(vvmax)/float(tailsupport))
    logger.info('outer regret {}  (due to samples: {} due to tail: {}'.format(racc,tmp,racc-tmp))
    #regret lower bound
    rlow=0.
    for i,y in enumerate(sp.linspace(Yout[0]-nstd*sp.sqrt(vvmax),mvmax,tailsupport)):
        rlow+= gpbo.core.GPdc.EI(-y,-m,v)[0,0]*sp.stats.norm.pdf(y,mvmax,sp.sqrt(vvmax))*((mvmax-Yout[0]+nstd*sp.sqrt(vvmax))/float(tailsupport))

    #regret from samples
    rsam=0.
    for i in xrange(Q.shape[0]):
        rsam+=max(0.,Q[i,1]-Q[i,2])
    rsam/=Q.shape[0]

    #local regret from incumbent from samples
    rloc=0.
    for i in xrange(Q.shape[0]):
        rloc+=max(0.,Q[i,3]-Q[i,1])
    rloc/=Q.shape[0]

    #set switch to local if condition achieved
    rval=0
    if racc<para['regretswitch']:
            rval=1
            persist['flip']=True
            optstate.startlocal=xmin
    if gpbo.core.debugoptions['adaptive']:
        if d==2:
            gpbo.core.optutils.plotaslogrtheta(W[:,0],W[:,1],xmin[0],xmin[1],ax[1,1],'b.')
            ax[0,2].plot(W[:,0],W[:,1],'b.')
            #plot mindraws
            R, Y, A = ESutils.draw_min_xypairgrad(G, W, 1500, xmin)
            gpbo.core.optutils.plotaslogrtheta(R[:,0],R[:,1],xmin[0],xmin[1],ax[1,1],'r.')
            ax[0,2].plot(R[:,0],R[:,1],'r.')
        ax[2,2].plot(Q[:,1],Q[:,2],'r.')
        ax[2,2].set_xlabel('inR')
        ax[2,2].set_ylabel('outR')
        ax[2,2].plot([ymin],[ymin],'go')

        ax[2,1].plot(Q[:,1],Q[:,3],'r.')
        ax[2,1].set_xlabel('inR')
        ax[2,1].set_ylabel('atArg')
        ax[2,1].plot([ymin],[ymin],'go')

        def pltcdf(Y,C,ax,col):
            return ax.plot(sp.hstack([[i,i] for i in Y]),sp.hstack([[i-C[0],i] for i in C]),col)

        pltcdf(Yin,Cin,ax[2,0],'b')
        rin = sp.linspace(Yin[0],Yin[-1],150)
        ax[2,0].plot(rin, map(lambda x:sp.stats.norm.cdf(x,*normin),rin),'k')


        pltcdf(Yat,Cat,ax[2,0],'g')
        rat = sp.linspace(Yat[0],Yat[-1],150)
        ax[2,0].plot(rat, map(lambda x:sp.stats.norm.cdf(x,*normat),rat),'k')
        ax[2,0].set_yscale('logit')

        pltcdf(Yout,Cout,ax[1,0],'r')
        ax[1,0].set_yscale('logit')

        mxo=Yout[-1]
        mno=Yout[0]
        ro = sp.linspace(min(mno-0.05*(mxo-mno),ymin),mxo+0.05*(mxo-mno),200)

        ax[1,0].plot(ro,map(lambda x:sp.stats.norm.cdf(x,loc=mvmax,scale=sp.sqrt(vvmax)),ro),'g')

        ax[1,0].plot(ro,map(lambda x:sp.stats.norm.cdf(x,loc=mu,scale=sp.sqrt(vvmax)),ro),'b',alpha=0.5)
        ax[1,0].plot(ro,map(binormdist.cdf,ro),'k')
        ax[1,0].plot(ro,map(binormdistF.cdf,ro),'purple')
        ax[1,2].text(0,0.2,'regretest{}'.format(racc))
        ax[1,2].text(0,0.15,'regretsam{}'.format(rsam))
        ax[1,2].text(0,0.3,'localrsam{}'.format(rloc))
        ax[1,2].text(0,0.1,'regretlow {} '.format(rlow))
        ax[1,2].text(0,0.4,'localrest {} '.format(lrest))
        persist['Rexists'].append(optstate.n)
        persist['sampleregret'].append(rsam)
        persist['expectedregret'].append(racc)
        persist['localrsam'].append(rloc)
        persist['regretlower'].append(rlow)
        persist['localrest'].append(lrest)
        ax[0,3].plot(persist['Rexists'],persist['localrest'],'k')
        ax[0,3].plot(persist['Rexists'],persist['sampleregret'],'b')
        ax[0,3].plot(persist['Rexists'],persist['expectedregret'],'g')
        ax[0,3].plot(persist['Rexists'],persist['localrsam'],'r')
        #ax[0,3].plot(persist['Rexists'],persist['regretlower'],'purple')
        ax[0,3].set_yscale('log')
        ax[2,3].set_yscale('log')
        try:
            from gpbo.core import debugpath
            fname = 'lotsofplots' + time.strftime('%d_%m_%y_%H:%M:%S') + '.png'
            print('saving as {}'.format(fname))
            fig.savefig(os.path.join(debugpath, fname))
        except BaseException as e:
            logger.error(str(e))
        fig.clf()
        plt.close(fig)
        del (fig)

    #if a cheat objective as available see how we would do on starting a local opt now
    if 'cheatf' in para.keys():
        try:
            C = sp.linalg.cholesky(H)
        except:
            logger.info('not +ve definite at posterior min')
            C=sp.linalg.cholesky(sp.eye(H.shape[0]))
        print('C {} \nxmin {}\nC.T.xmin{}'.format(C,xmin,C.T.dot(xmin)))
        def fn2(x):
            print(x,para['cheatf'](sp.linalg.solve(C.T,x),**{'s':0.,'d':[sp.NaN]})[0])
            return para['cheatf'](sp.linalg.solve(C.T,x),**{'s':0.,'d':[sp.NaN]})[0]
        R=minimize(fn2,C.T.dot(xmin),method='bfgs')
        logger.warn('cheat testopt result with precondition {}:\n{}'.format(H,R))

    return rval,persist,{'start':xmin,'H':H,'reuseH':[k.hyp for k in G.kf]}


def alternate(optstate,persist,**para):
    return optstate.n%2,None,dict()


def aftern(optstate,persist,**para):
    if optstate.n>para['n']:
        return 1,None,{'start':[0.,0.]}
    else:
        return 0,None,dict()


