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


    #pvecount=0
    #for i in xrange(2000):
    #    Hdraw=gpbo.core.optutils.Hvec2H(Hdist.rvs(),d)
    #    try:
    #        sp.linalg.cholesky(Hdraw)
    #        pvecount+=1
    #    except sp.linalg.LinAlgError:
    #        pass
    #pvecount/=2000.

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

    C = sp.linalg.cholesky(H)
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
        print( 'plotting choice...')
        fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(30, 30))
        #plot the current GP
        n = 50
        x_ = sp.linspace(-1, 1, n)
        y_ = sp.linspace(-1, 1, n)
        z_ = sp.empty([n, n])
        s_ = sp.empty([n, n])
        for i in range(n):
            for j in range(n):
                m_, v_ = G.infer_diag_post(sp.array([y_[j], x_[i]]), [[sp.NaN]])
                z_[i, j] = m_[0, 0]
                s_[i, j] = sp.sqrt(v_[0, 0])
        print('0')
        CS = ax[0, 0].contour(x_, y_, z_, 20)
        ax[0, 0].clabel(CS, inline=1, fontsize=10)
        print('1')
        CS = ax[1, 0].contour(x_, y_, s_, 20)
        ax[0, 0].axis([-1., 1., -1., 1.])
        ax[1, 0].clabel(CS, inline=1, fontsize=10)
        print('2')
        ax[0, 0].plot(xmin[0], xmin[1], 'ro')


        #plot support points and draws
        print('3')
#        for i in tqdm.tqdm(xrange(para['support'])):
#            symbol = ['r.','g.','bx'][Z_[i]]
#            ax[0,1].plot(W[i,0],W[i,1],symbol)
#            if ER[0,i]>0.:
#                ax[2,2].plot(W[i,0],W[i,1],'b.')
#            else:
#                ax[2,2].plot(W[i,0],W[i,1],'r.')
#        ax[0, 1].axis([-1, 1, -1, 1])

#        print('5')
#        for i in tqdm.tqdm(xrange(nd)):
#            symbol = ['r.', 'g.'][Y_[i]]
#            ax[1, 1].plot(R[i, 0], R[i, 1], symbol)
#        for i in xrange(nd):
#            symbol = ['r.', 'g.'][Y_r[i]]
#            ax[1, 2].plot(R[i, 0], R[i, 1], symbol)

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
        ax[2,0].plot(persist['ERegret'],'k') #full regret
        ax[2, 0].plot(persist['LRegret'], 'b') #local regret
        ax[2, 0].plot(persist['GRegret'], 'm') #nonlocal regret


        #plot support based regret bounds
        ax[2, 0].plot(persist['EBound'], 'k--.')  # full regret
        ax[2, 0].plot(persist['LBound'], 'b--.')  # local regret
        ax[2, 0].plot(persist['GBound'], 'm--.')  # nonlocal regret


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
        ax[2, 0].plot(xaxis,Lp, 'b:')  # local regret
        ax[2, 0].plot(xaxis,Gp, 'm:')


        #true regret if available:
        #try:
        Rtrue = para['cheatf'](xmin,**{'s':0.,'d':[sp.NaN]})[0]-para['cheatymin']
        persist['RTrue'].append(Rtrue)
        ax[2,1].plot(persist['RTrue'],'g')
        ax[2,1].plot(persist['ERegret'],'k') #prob not in local

        try:
            ax[2,0].set_yscale('log')
        except:
            pass
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
        print(grf(i-1)+lsr,grf(i-1)+lsr<tol,tol,Sc[1,i],Sc[1,i]<brm)
        Sv[1,i]=grf(i-1)+lsr+Sc[1,i] if Sc[1,i]<brm and grf(i-1)+lsr<tol else sp.Inf
        if Sv[1,i+1]>=Sv[0,i+1]:
            Sv[0,i]=Sv[0,i+1]
            Ac[i]=0
        else:
            Sv[0,i]=Sv[1,i+1]
            Ac[i]=1
            switchat=i
    return Ac[0],switchat