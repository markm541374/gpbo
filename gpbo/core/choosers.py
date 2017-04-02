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

def introspection(optstate,persist,**para):
    logging.info('\n--------------------------------\nIntrospection\n')
    if persist == None:
        persist = defaultdict(list)
    if optstate.n < para['onlyafter']:
        return 0, persist, dict()
    if persist['flip']:
        return 1, persist, dict()
    rval=0
    fig, ax = plt.subplots(nrows=3, ncols=4, figsize=(85, 85))
    d = len(para['lb'])
    lb = para['lb']
    ub = para['ub']
    #build a GP with slice-samples hypers
    x = sp.vstack(optstate.x)
    y = sp.vstack(optstate.y)
    s = sp.vstack([e['s'] for e in optstate.ev])
    dx = [e['d'] for e in optstate.ev]
    G = PES.makeG(x, y, s, dx, para['kindex'], para['mprior'], para['sprior'], para['nhyp'])


    # plot the current GP
    gpbo.core.optutils.gpplot(ax[0,0],ax[0,1],G,para['lb'],para['ub'],ns=60)
    ax[0,0].set_title('GP_post_mean')
    ax[0,1].set_title('GP_post_var')

    #find the pred. minimum
    #[dxmin,dymin,dierror] = direct(lambda x:G.infer_m_post(x,[[sp.NaN]])[0,0],para['lb'],para['ub'],user_data=[], algmethod=1, maxf=8000, logfilename='/dev/null')
    #xmin,ymin,ierror= gpbo.core.optutils.boundedlocal(lambda x:G.infer_m_post(x,[[sp.NaN]])[0,0],para['lb'],para['ub'],dxmin,**{'gtol':0.00001,'maxfun':250})
    dpara = {'user_data':[],
             'algmethod':1,
             'maxf':8000,
             'logfilename':'/dev/null'}
    lpara = {'gtol':0.00001,
             'maxfun':300}

    xmin,ymin,ierror = gpbo.core.optutils.twopartopt(lambda x:G.infer_m_post(x,[[sp.NaN]])[0,0],para['lb'],para['ub'],dpara,lpara)
    xvmax,vmax,ierror = gpbo.core.optutils.twopartopt(lambda x:-G.infer_diag_post(x,[[sp.NaN]])[1][0,0],para['lb'],para['ub'],dpara,lpara)
    mvmax,vvmax = [j[0,0] for j in G.infer_diag_post(xvmax,[[sp.NaN]])]
    print('\n----------------------VMAX {}\n {} {}'.format(xvmax,mvmax,vvmax))
    ax[0, 0].plot(xmin[0], xmin[1], 'ro')

    #get hessian/grad posterior
    #local probmin elipse at post min
    GH = gpbo.core.optutils.gpGH(G,xmin)
    Gr,cG,H,Hvec,varHvec,_,_ = GH
    Hdist = sp.stats.multivariate_normal(Hvec.flatten(),varHvec)


    #plot some draws from H
    for i in xrange(20):
        Gm,Gv,Hd = gpbo.core.drawconditionH(*GH)
        try:
            sp.linalg.cholesky(Hd)
            gpbo.core.optutils.plotprobstatellipse(Gv,Hd,xmin,ax[1,1],logr=True)
        except sp.linalg.LinAlgError:
            pass

    #ax[1,0].text(xmin[0],xmin[1],'{}'.format(pvecount,))
    ax[1,1].set_title('local pmin ellipse')

    #step out checking for +ve definite
    rmax=0
    done=False
    Rad = sp.logspace(-4,0,200)
    PP = sp.empty(Rad.size)
    for j,r in tqdm.tqdm(enumerate(Rad)):
        theta = sp.random.uniform(0,2*sp.pi)
        x0 = xmin[0]+sp.sin(theta)*r
        x1 = xmin[1]+sp.cos(theta)*r
        p = gpbo.core.optutils.probgppve(G,sp.array([x0,x1]))
        PP[j]=p
        if p<1-1e-5 and not done:
            if j>0:
                rmax=Rad[j-1]
            else:
                rmax=0
            done=True
    circ = sp.empty([2, 100])
    for i in xrange(100):
        theta = 2. * sp.pi * i / 99.
        circ[:, i]=[xmin[0]+rmax*sp.sin(theta),xmin[1]+rmax*sp.cos(theta)]
    if rmax>0:
        ax[1,1].plot([sp.log10(rmax)]*2,[0.,2*sp.pi],'purple')
    ax[1,3].set_title('pve definite region')
    ax[1,3].semilogx(Rad,PP)

    #draw support points
    W = sp.vstack(ESutils.draw_support(G, lb, ub, 2500, ESutils.SUPPORT_LAPAPROT, para=20))
    #add some points loguniform around -3 to -1 to see if it makes a difference
    Waug = sp.empty([500,2])
    for i in xrange(500):
        r = sp.random.uniform(-4,0)
        theta = sp.random.uniform(0,sp.pi*2)
        Waug[i,0] = min(1.,max(-1.,xmin[0]+sp.sin(theta)*10**r))
        Waug[i,1] = min(1.,max(-1.,xmin[1]+sp.cos(theta)*10**r))
    W = sp.vstack([W,Waug])
    nd = 1500

    #draw mins and value of g at xmin as pair
    R, Y, A = ESutils.draw_min_xypairgrad(G, W, nd, xmin)
    #plot support
    gpbo.core.optutils.plotaslogrtheta(W[:,0],W[:,1],xmin[0],xmin[1],ax[1,1],'b.')
    ax[0,2].plot(W[:,0],W[:,1],'b.')
    #plot mindraws
    gpbo.core.optutils.plotaslogrtheta(R[:,0],R[:,1],xmin[0],xmin[1],ax[1,1],'r.')
    ax[0,2].plot(R[:,0],R[:,1],'r.')





    if rmax==0:
        ax[1,2].text(0,0,'no pregion')
    else:
        #targeteps = 1e-5
        #targetn = int((1 - 2 * targeteps) / float(targeteps))
        #batchsize = 5000
        #nbatches = 1 + targetn / batchsize
        #for i in tqdm.tqdm(xrange(nbatches)):
        #    R, Y, A = ESutils.draw_min_xypairgrad(G, W, batchsize, xmin)
        #    count = 0
        #    for j in range(batchsize):
        #        if sp.sqrt((R[j, 0] - xmin[0]) ** 2 + (R[j, 1] - xmin[1]) ** 2) > rmax:
        #            count += 1
        #    if count > 0:
        #        break
        #nout = count
        #nin = i * batchsize + batchsize - count
        #pinregion = (nin+1) / float(nin + nout+2)
        #ax[1,2].text(0,0.5,'pregionsize {}\np(x* !in R)={} (from {} mindraws) '.format(rmax,1.-pinregion,nin+nout))

        Q = gpbo.core.optutils.drawpartitionmin(G,W,xmin,rmax,8000)
        nin = sum(Q[:,4])
        #ax[1,2].text(0,0.75,'{}'.format(1-nin/2000.))
        ax[2,2].plot(Q[:,1],Q[:,2],'r.')
        ax[2,2].set_xlabel('inR')
        ax[2,2].set_ylabel('outR')
        ax[2,2].plot([ymin],[ymin],'go')

        ax[2,1].plot(Q[:,1],Q[:,3],'r.')
        ax[2,1].set_xlabel('inR')
        ax[2,1].set_ylabel('atArg')
        ax[2,1].plot([ymin],[ymin],'go')

        #pcurves from Q
        def data2cdf(X):
            n = X.size
            C = sp.linspace(1./n,1,n)
            XC = sorted(X)
            #mn = sp.mean(XC)
            #var = sp.var(XC)
            #axc.plot(XC,map(lambda x:sp.stats.norm.cdf(x,loc=mn,scale=sp.sqrt(var)),XC),'k')
            #axc.plot(XC[1:],C[1:],col+'.-')
            return XC,C
        def pltcdf(Y,C,ax,col):
            return ax.plot(sp.hstack([[i,i] for i in Y]),sp.hstack([[i-C[0],i] for i in C]),col)



        Yin,Cin = data2cdf(Q[:,1])
        #ax[2,0].plot(Yin,Cin,'b.')
        pltcdf(Yin,Cin,ax[2,0],'b')
        normin = sp.stats.norm.fit(Yin)
        rin = sp.linspace(Yin[0],Yin[-1],150)
        ax[2,0].plot(rin, map(lambda x:sp.stats.norm.cdf(x,*normin),rin),'k')


        Yat,Cat = data2cdf(Q[:,3])
        #ax[2,0].plot(Yat,Cat,'g.')
        pltcdf(Yat,Cat,ax[2,0],'g')
        normat = sp.stats.norm.fit(Yat)
        rat = sp.linspace(Yat[0],Yat[-1],150)
        ax[2,0].plot(rat, map(lambda x:sp.stats.norm.cdf(x,*normat),rat),'k')

        ax[2,0].set_yscale('logit')

        Yout,Cout = data2cdf(Q[:,2])
        #ax[1,0].plot(Yout,Cout,'r.-')
        pltcdf(Yout,Cout,ax[1,0],'r')

        ax[1,0].set_yscale('logit')

        mxo=Yout[-1]
        mno=Yout[0]
        ro = sp.linspace(min(mno-0.05*(mxo-mno),ymin),mxo+0.05*(mxo-mno),200)

        ax[1,0].plot(ro,map(lambda x:sp.stats.norm.cdf(x,loc=mvmax,scale=sp.sqrt(vvmax)),ro),'g')

      #  ax[1,0].axis([mno-0.05*(mxo-mno),mxo+0.05*(mxo-mno),0.05/Q.shape[0],1-0.1/Q.shape[0]])

        mi,vi = sp.stats.norm.fit(Q[:,1])
        pin = 1.- sp.stats.norm.cdf(0.,loc=mi-mvmax,scale=sp.sqrt(vvmax+vi**2))
        ax[1,2].text(0,0.9,'pign {}'.format(pin))

        ymin=Yout[0]
        cdfymin=Cout[0]
        mu = ymin-sp.sqrt(vvmax*2)*sp.special.erfinv(2*cdfymin-1.)
        ax[1,0].plot(ro,map(lambda x:sp.stats.norm.cdf(x,loc=mu,scale=sp.sqrt(vvmax)),ro),'b',alpha=0.5)
        #ax[1,0].plot(ymin,cdfymin,'ro')
        print('\n-----------------------\ny{} c{}\ninfered c {}'.format(ymin,cdfymin,sp.stats.norm.cdf(ymin,loc=mu,scale=sp.sqrt(vvmax))))
        pin2 = 1.-sp.stats.norm.cdf(0.,loc=mi-mu,scale=sp.sqrt(vvmax+vi**2))

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
        #ax[0,3].plot(ro,map(splicecdf,ro),'r')
        #ax[0,3].set_yscale('log')
        #prob of less than in
        acc=0
        for i in xrange(150):
            yindraw = sp.stats.norm.rvs(*normin)
            punderydraw = splicecdf(yindraw)
            acc+=punderydraw
        acc/=150.
        pestupper = acc
        #try:

        #expected regrets
        racc = 0.
        m,v=normin
        n=len(Cout)
        for i in xrange(1,n):
            racc+= gpbo.core.GPdc.EI(-Yout[i],-m,sp.sqrt(v))[0,0]/float(n)
        #print('racc{}'.format(racc))
        for i,y in enumerate(sp.linspace(Yout[0]-3*sp.sqrt(vvmax),Yout[0],100)):
            racc+= gpbo.core.GPdc.EI(-y,-m,sp.sqrt(v))[0,0]*sp.stats.norm.pdf(y,mu,sp.sqrt(vvmax))*(3*sp.sqrt(vvmax)/100.)
            #print('y{} EI{} p{} r{}'.format(y, gpbo.core.GPdc.EI(-y,-m,sp.sqrt(v)),sp.stats.norm.pdf(y,mu,sp.sqrt(vvmax)),gpbo.core.GPdc.EI(-y,-m,sp.sqrt(v))*(3*sp.sqrt(vvmax)/100.)*sp.stats.norm.pdf(y,mu,sp.sqrt(vvmax))))

        rlow=0.
        for i,y in enumerate(sp.linspace(Yout[0]-3*sp.sqrt(vvmax),mvmax,200)):
            rlow+= gpbo.core.GPdc.EI(-y,-m,sp.sqrt(v))[0,0]*sp.stats.norm.pdf(y,mvmax,sp.sqrt(vvmax))*((mvmax-Yout[0]+3*sp.sqrt(vvmax))/200.)

        rsam=0.
        for i in xrange(Q.shape[0]):
            rsam+=max(0.,Q[i,1]-Q[i,2])
        rsam/=Q.shape[0]


        rloc=0.
        for i in xrange(Q.shape[0]):
            rloc+=max(0.,Q[i,3]-Q[i,1])
        rloc/=Q.shape[0]

        if racc<para['regretswitch']:
            rval=1
            persist['flip']=True
            optstate.startlocal=xmin
        ax[1,2].text(0,0.2,'regretest{}'.format(racc))
        ax[1,2].text(0,0.15,'regretsam{}'.format(rsam))
        ax[1,2].text(0,0.3,'localrsam{}'.format(rloc))
        ax[1,2].text(0,0.8,'pestupper {}'.format(pestupper))
        ax[1,2].text(0,0.1,'regretlow {} '.format(rlow))
        persist['Rexists'].append(optstate.n)
        persist['predictedpinR'].append(pin)
        persist['higherpredictedpinR'].append(pestupper)
        #persist['sampledpinR'].append(1.-pinregion)
        persist['sampleregret'].append(rsam)
        persist['expectedregret'].append(racc)
        persist['localrsam'].append(rloc)
        persist['regretlower'].append(rlow)
        ax[0,3].plot(persist['Rexists'],persist['sampleregret'],'b')
        ax[0,3].plot(persist['Rexists'],persist['expectedregret'],'g')
        ax[0,3].plot(persist['Rexists'],persist['localrsam'],'r')
        ax[0,3].plot(persist['Rexists'],persist['regretlower'],'purple')
        ax[0,3].set_yscale('log')
        ax[2,3].plot(persist['Rexists'],persist['predictedpinR'],'r')
        ax[2,3].plot(persist['Rexists'],persist['higherpredictedpinR'],'g')
        #ax[2,3].plot(persist['Rexists'],persist['sampledpinR'],'b')
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
    def fn(x):
        return para['cheatf'](x,**{'s':0.,'d':[sp.NaN]})[0]
    R=minimize(fn,xmin,method='bfgs')
    logger.warn('cheat testopt result from {}:\n{}'.format(xmin,R))

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
    return rval,persist,{'start':xmin,'H':H}


def alternate(optstate,persist,**para):
    return optstate.n%2,None,dict()


def aftern(optstate,persist,**para):
    if optstate.n>para['n']:
        return 1,None,{'start':[0.,0.]}
    else:
        return 0,None,dict()


