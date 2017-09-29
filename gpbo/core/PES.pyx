# Classes to run PES. Classes for constant or variable obs noise and augmented space

#make a single posterior gp form data and take draws on this
from __future__ import print_function
xrange=range
import sys
from gpbo.core import ESutils
from gpbo.core import GPdc
from gpbo.core import eprop
import scipy as sp
from scipy import stats as sps
#from gpbo.core.optutils import silentdirect as direct

try:
    from matplotlib import pyplot as plt
    plots=True
except ImportError:
    plots=False
    plt=None

import logging
from libc.math cimport log10, log, exp, M_E

logger = logging.getLogger(__name__)
import gpbo



import os
import time

def makeG(X,Y,S,D,kindex,mprior,sprior,nh,chains=1,prior='lognorm'):
    #draw hyps based on plk
    #print "RRRRRRRRRRRRRR"+str([X,Y,S,D,kindex,mprior,sprior,nh])
    H = ESutils.drawhyp_plk(X,Y,S,D,kindex,mprior,sprior,nh,chains=chains,prior=prior)
    
    G = GPdc.GPcore(X, Y, S, D, [GPdc.kernel(kindex, X.shape[1], i) for i in H])
    
    return G

def drawmins(G,n,lb,ub,SUPPORT=300,mode = [ESutils.SUPPORT_SLICELCB],SLICELCB_PARA=1.):
    #draw support points
    W = sp.vstack([ESutils.draw_support(G, lb,ub,SUPPORT/len(mode),m, para = SLICELCB_PARA) for m in mode])

    R = ESutils.draw_min(G,W,n)
    #draw in samples on the support
    print( "drawing mins from support")

    from gpbo.core import debugoutput
    if debugoutput['support'] and plots:
        fig, ax = plt.subplots(1)
        ax.plot(W[:,0],W[:,1],'bx')
        ax.plot(R[:,0],R[:,1],'r.')
        from gpbo.core import debugpath
        fig.savefig(os.path.join(debugpath,'support'+time.strftime('%d_%m_%y_%H:%M:%S')+'.png'))
        fig.clf()
        plt.close(fig)
        del(fig)
    #plt.show()
    return R

def drawmins_inplane(G,n,lb,ub,axis,value, SUPPORT=300, mode=ESutils.SUPPORT_SLICELCB, SLICELCB_PARA=1.):
    W = sp.vstack([ESutils.draw_support_inplane(G, lb,ub,SUPPORT/len(mode),m, axis,value,para = SLICELCB_PARA) for m in mode])
    #draw in samples on the support
    R = ESutils.draw_min(G,W,n)
    return R

OFFHESSZERO=0
OFFHESSINFER=1

def addmins(G,X,Y,S,D,xmin,mode=OFFHESSZERO, GRADNOISE=1e-9,EP_SOFTNESS=1e-9,int EPROP_LOOPS=20,dropedge=True):
    cdef int i
    cdef int dim=X.shape[1]
    dropdims=[]
    for i in range(dim):
        if xmin[i]>0.99 or xmin[i]<-0.99:
            dropdims.append(i)
    #grad elements are zero
    Xg = sp.vstack([xmin]*dim)
    Yg = sp.zeros([dim,1])
    Sg = sp.ones([dim,1])*GRADNOISE
    Dg = [[i] for i in range(dim)]
    if dropedge:
        for i in dropdims:
            Sg[i,0]=1e9

    #offdiag hessian elements
    nh = ((dim-1)*dim)/2
    Xh = sp.vstack([sp.empty([0,dim])]+[xmin]*nh)
    class MJMError(Exception):
        pass
    if mode==OFFHESSZERO:
        Yh = sp.zeros([nh,1])
        Sh = sp.ones([nh,1])*GRADNOISE


    elif mode==OFFHESSINFER:
        raise MJMError("addmins with mode offhessinfer not implemented yet")
    else:
        raise MJMError("invalid mode in addmins")
    Dh=[]
    for i in xrange(dim):
        for j in xrange(i):
            Dh.append([i,j])
            if (i in dropdims or j in dropdims) and dropedge:
               Sh[len(Dh)-1]=1e9
    #diag hessian and min
    Xd = sp.vstack([xmin]*(dim+1))
    Dd = [[sp.NaN]]+[[i,i] for i in xrange(dim)]
    [m,V] = G.infer_full_post(Xd,Dd)
    for i in xrange(dim):
        if V[i,i]<0:
            class MJMError(Exception):
                pass
            print( [m,V])
            raise MJMError('negative on diagonal')
        
    yminarg = sp.argmin(Y)
    Y_ = sp.array([Y[yminarg,0]]+[0.]*dim)
    Z = sp.array([-1]+[1.]*dim)
    F = sp.array([S[yminarg,0]]+[EP_SOFTNESS]*dim)
    [Yd,Stmp] = eprop.expectation_prop(m,V,Y_,Z,F,EPROP_LOOPS)
    Sd = sp.diagonal(Stmp).flatten()
    Sd.resize([dim+1,1])
    Yd.resize([dim+1,1])
    if dropedge:
        for i in dropdims:
            Sd[i,0]=1e9

    #concat the obs
    Xo = sp.vstack([X,Xg,Xd,Xh])
    Yo = sp.vstack([Y,Yg,Yd,Yh])
    So = sp.vstack([S,Sg,Sd,Sh])
    Do = D+Dg+Dd+Dh
    
    return [Xo,Yo,So,Do]

NOMIN=0
def addmins_inplane(G,X,Y,S,D,xmin,axis,value,mode=OFFHESSZERO, GRADNOISE=1e-9,EP_SOFTNESS=1e-9,EPROP_LOOPS=20,MINPOLICY=NOMIN,dropedge=True):
    cdef int i,j
    cdef int dim=X.shape[1]

    #dropdims=[]
    #for i in range(dim):
    #    if (xmin[i]>0.99 or xmin[i]<-0.99) and i!=axis:
    #        dropdims.append(i)
    #grad elements are zero
    Xg = sp.vstack([xmin]*(dim-1))
    Yg = sp.zeros([dim-1,1])
    Sg = sp.ones([dim-1,1])*GRADNOISE
    Dg = [[i] for i in range(dim) if i!=axis]

    #if dropedge:
    #    for i in dropdims:
    #        Sg[i,0]=1e9
    #offdiag hessian elements
    nh = ((dim-1)*(dim-2))/2
    Xh = sp.vstack([sp.empty([0,dim])]+[xmin]*nh)
    Dh=[]
    for i in range(dim):
        for j in range(i):
            if i!=axis and j!=axis:
                Dh.append([i,j])
    class MJMError(Exception):
        pass
    if mode==OFFHESSZERO:
        Yh = sp.zeros([nh,1])
        Sh = sp.ones([nh,1])*GRADNOISE
    elif mode==OFFHESSINFER:
        raise MJMError("addmins with mode offhessinfer not implemented yet")
    else:
        raise MJMError("invalid mode in addmins")
    #diag hessian and min
    if MINPOLICY==NOMIN:
        Xd = sp.vstack([xmin]*(dim-1))
        Dd = [[i,i] for i in xrange(dim) if i !=axis]
        [m,V] = G.infer_full_post(Xd,Dd)
        #yminarg = sp.argmin(Y)
    
        Y_ = sp.array([0.]*dim)
        Z = sp.array([1.]*dim)
        F = sp.array([EP_SOFTNESS]*dim)

        [Yd,Stmp] = eprop.expectation_prop(m,V,Y_,Z,F,EPROP_LOOPS)
        Sd = sp.diagonal(Stmp).flatten()
        Sd.resize([dim-1,1])
        Yd.resize([dim-1,1])
    else:
        print( "this isn't a valid approach!!!!!!!!!")
        Xd = sp.vstack([xmin]*(dim-1+1))
        Dd = [[sp.NaN]]+[[i,i] for i in xrange(dim) if i !=axis]
        [m,V] = G.infer_full_post(Xd,Dd)
        yminarg = sp.argmin(Y)
    
        Y_ = sp.array([Y[yminarg,0]]+[0.]*dim)
        Z = sp.array([-1]+[1.]*dim)
        F = sp.array([S[yminarg,0]]+[EP_SOFTNESS]*dim)

        [Yd,Stmp] = eprop.expectation_prop(m,V,Y_,Z,F,EPROP_LOOPS)
        Sd = sp.diagonal(Stmp).flatten()
        Sd.resize([dim+1-1,1])
        Yd.resize([dim+1-1,1])
    #concat the obs
    Xo = sp.vstack([X,Xg,Xd,Xh])
    Yo = sp.vstack([Y,Yg,Yd,Yh])
    So = sp.vstack([S,Sg,Sd,Sh])
    Do = D+Dg+Dd+Dh
    
    return [Xo,Yo,So,Do]

def PESgain(g0,G1,Z,X,D,s):
    
    H = sp.zeros(len(D))
    [m0,v0] = g0.infer_diag_post(X,D)
    
    #print X.shape
    for j in xrange(X.shape[0]):
        H[j]-= len(G1)*0.5*log(2*sp.pi*M_E*(v0[0,j]+s[j]))
        for i,g1 in enumerate(G1):
            
            Xi = sp.vstack([X[j,:],Z[i,:]])
            Di = [D[j]]+[[sp.NaN]]
            
            [mi,Vi] = g1.infer_full_post(Xi,Di)
            #print [mi,Vi]
            v1 = Vadj(mi,Vi)
        
            h1 = 0.5*log(2*sp.pi*M_E*(v1+s[j]))
            H[j]+=h1
    
    return -H/float(len(G1))

def Vadj(m,V):
    s = V[0,0]+V[1,1]-2*V[0,1]
    mu = m[0,1]-m[0,0]
    alpha = -mu/sp.sqrt(s) #sign difference for minimizing
    try:
        beta = exp(sps.norm.logpdf(alpha) - sps.norm.logcdf(alpha))
    except:
        logger.warn('Vadj: complex in PESgain, using beta=0\n alpha: {}\nmu: {}\ns: {}\nm: {}\nV: {}\n'.format(alpha,mu,s,m,V))
        beta=0.
    #print [s,mu,alpha,beta]
    vadj = V[0,0]-beta*(beta+alpha)*(1./s)*(V[0, 0]-V[0, 1])**2
    return vadj


#basic PES class if search_pes is used. variable noise if search_acq is used
class PES:
    def __init__(self,X,Y,S,D,lb,ub,kindex,mprior,sprior,DH_SAMPLES=8,DM_SAMPLES=8, DM_SUPPORT=400,DM_SLICELCBPARA=1.,mode=ESutils.SUPPORT_SLICELCB,noS=False,DM_DROP=True,preselectH=False):
        print( "PES init:")
        self.lb=lb
        self.ub=ub
        self.noS=noS
        if noS:
            S=sp.zeros(S.shape)
        if not preselectH:
            self.G = makeG(X,Y,S,D,kindex,mprior,sprior,DH_SAMPLES)
        else:
            logger.info('reusing preselected hyperparameters')
            self.G =  GPdc.GPcore(X,Y,S,D, [GPdc.kernel(kindex, X.shape[1], h) for h in preselectH])
        HS = sp.vstack([k.hyp for k in self.G.kf])
        self.Z = drawmins(self.G,DM_SAMPLES,lb,ub,SUPPORT=DM_SUPPORT,SLICELCB_PARA=DM_SLICELCBPARA,mode=mode)
        #print "mindraws: "+str(self.Z)
        self.Ga = [GPdc.GPcore(*addmins(self.G, X, Y, S, D, self.Z[i, :],dropedge=DM_DROP) + [self.G.kf]) for i in xrange(DM_SAMPLES)]
    def __del__(self):
        try:
            self.G.__del__()
        except:
            pass
        try:
            self.Ga.__del__()
        except:
            pass
        return    
    def query_pes(self,Xq,Sq,Dq):
        
        a = PESgain(self.G,self.Ga,self.Z,Xq,Dq,Sq)
        return a
    
    def query_acq(self,Xq,Sq,Dq,costfn):
        a = PESgain(self.G,self.Ga,self.Z,Xq,Dq,Sq)
        for i in xrange(Xq.shape[0]):
            a[i] = a[i]/costfn(Xq[i,:].flatten(),Sq[i,:].flatten())
        return a



    def search_vmax(self,s,spara,dv=[[sp.NaN]],):
        self.stmp = s
        def wrap(Q):
            x = sp.array([Q])
            v = self.G.infer_diag(x,[[sp.NaN]])[1][0,0]
            return -v
        #logger.critical('aaaa')
        xmin,ymin,ierror = gpbo.core.optutils.twopartopt(wrap,self.lb,self.ub,spara['dpara'],spara['lpara'])
        #[xmin, ymin, ierror] = direct(directwrap,self.lb,self.ub,user_data=[], algmethod=1, maxf=maxf, logfilename='/dev/null')

        return [xmin,-ymin,ierror]

    def search_pes(self,s,spara,dv=[[sp.NaN]],):
        self.stmp = s
        def wrap(Q):
            x = sp.array([Q])
            if self.noS:
                alls = [k(x,x,dv,dv,gets=True)[1] for k in self.G.kf]
                s = exp(sp.mean(log(alls)))
            else:
                s= self.stmp
            acq = PESgain(self.G,self.Ga,self.Z,x,dv,[s])
            R = -acq
            return R
        #logger.critical('aaaa')
        xmin,ymin,ierror = gpbo.core.optutils.twopartopt(wrap,self.lb,self.ub,spara['dpara'],spara['lpara'])
        #[xmin, ymin, ierror] = direct(directwrap,self.lb,self.ub,user_data=[], algmethod=1, maxf=maxf, logfilename='/dev/null')

        return [xmin,ymin,ierror]
    
    def search_acq(self,cfn,logsl,logsu,spara,dv=[[sp.NaN]],over=0.):
        print('over {}'.format(over))
        def wrap(Q):
            x = sp.array([Q[:-1]])
            s = 10**Q[-1]
            acq = PESgain(self.G,self.Ga,self.Z,x,dv,[s])
            try:
                R = -acq/(cfn(x,**{'s':s})+over)
            except TypeError:
                R = -acq/(cfn(x,s)+over)
            return R

        #[xmin, ymin, ierror] = direct(directwrap,sp.hstack([self.lb,logsl]),sp.hstack([self.ub,logsu]),user_data=[], algmethod=1, maxf=maxf, logfilename='/dev/null')
        xmin,ymin,ierror = gpbo.core.optutils.twopartopt(wrap,sp.hstack([self.lb,logsl]),sp.hstack([self.ub,logsu]),spara['dpara'],spara['lpara'])

        if gpbo.core.debugoutput['acqfn1d']:
            print( 'plotting acq1d...')
            import time
            from matplotlib import pyplot as plt
            f,a=plt.subplots(3,figsize=[8,10],sharex=True)
            nn = 100
            cost = sp.empty(nn)
            infgain = sp.empty(nn)
            acqfn = sp.empty(nn)
            srange = sp.linspace(logsl,logsu+4,nn)
            for i in xrange(nn):
                s=10**srange[i]
                cost[i] = cfn(xmin,**{'s':s})+over
                infgain[i] = PESgain(self.G,self.Ga,self.Z,sp.array([xmin[:-1]]),dv,[s])
                acqfn[i] = infgain[i]/cost[i]

            m,v = self.G.infer_diag_post(sp.array([xmin[:-1]]),[[sp.NaN]])
            a[1].plot([sp.log10(v[0,0]),sp.log10(v[0,0])],[0.01,1.],'r')
            a[0].plot(srange,cost)
            a[1].plot(srange,infgain)
            a[2].plot(srange,acqfn)
            a[0].set_yscale('log')
            a[1].set_yscale('log')
            a[2].set_yscale('log')
            f.savefig(os.path.join(gpbo.core.debugoutput['path'], 'aq1d' + time.strftime('%d_%m_%y_%H:%M:%S') + '.png'))
            plt.close(f)
            del(f)

        return [xmin,ymin,ierror]

#augmented space PES
class PES_inplane:
    def __init__(self,X,Y,S,D,lb,ub,kindex,mprior,sprior,axis,value,DH_SAMPLES=8,DM_SAMPLES=8, DM_SUPPORT=400,DM_SLICELCBPARA=1.,AM_POLICY=NOMIN,mode=ESutils.SUPPORT_SLICELCB,noS=False,DH_CHAINS=1):
        #print "PES init:"
        self.lb=lb
        self.ub=ub
        self.noS=noS
        self.X=X
        if noS:
            S=sp.zeros(S.shape)
        #print [X.shape,Y.shape,mprior,sprior]
        self.G = makeG(X,Y,S,D,kindex,mprior,sprior,DH_SAMPLES,chains=DH_CHAINS)
        HS = sp.vstack([k.hyp for k in self.G.kf])

        self.Z = drawmins_inplane(self.G,DM_SAMPLES,lb,ub,axis=axis,value=value,SUPPORT=DM_SUPPORT,SLICELCB_PARA=DM_SLICELCBPARA,mode=mode)
        #print "mindraws:\n"+str(self.Z)
        self.Ga = [GPdc.GPcore(*addmins_inplane(self.G, X, Y, S, D, self.Z[i, :], axis=axis, value=value, MINPOLICY=AM_POLICY) + [self.G.kf]) for i in xrange(DM_SAMPLES)]
        return
    
    def __del__(self):
        try:
            self.G.__del__()
        except:
            pass
        try:
            self.Ga.__del__()
        except:
            pass
        return
    def query_pes(self,Xq,Sq,Dq):
        a = PESgain(self.G,self.Ga,self.Z,Xq,Dq,Sq)
        return a
    
    def query_acq(self,Xq,Sq,Dq,costfn):                                                                                                            
        a = PESgain(self.G,self.Ga,self.Z,Xq,Dq,Sq)
        for i in xrange(Xq.shape[0]):
            a[i] = a[i]/costfn(Xq[i,:].flatten())
        return a
    
    def search_acq(self,cfn,sfn,spara,dv=[[sp.NaN]],over=0.):
        print( 'overhead={}'.format(over))
        def wrap(Q):
            sys.stdout.flush()
            x = sp.array([Q])
            if self.noS:
                alls = [k(x,x,dv,dv,gets=True)[1] for k in self.G.kf]
                s = exp(sp.mean(log(alls)))
            else:
                s = sfn(x)
            acq = PESgain(self.G,self.Ga,self.Z,x,dv,[s])
            try:
                R = -acq/(cfn(x[0,1:],**{'xa':x[0,0]})+over)
            except TypeError:
                R = -acq/(cfn(x,s)+over)
            return R
        #print self.ub

        xmin,ymin,ierror = gpbo.core.optutils.twopartopt(wrap,self.lb,self.ub,spara['dpara'],spara['lpara'])
        #[xmin, ymin, ierror] = direct(directwrap,self.lb,self.ub,user_data=[], algmethod=0, maxf=maxf, logfilename='/dev/null')
        
        
        if False and plots:

            import time
            D = self.lb.size
            ns=200
            f,a = plt.subplots(2*D)
            
            if self.noS:
                alls = [k(xmin,xmin,dv,dv,gets=True)[1] for k in self.G.kf]
                s = exp(sp.mean(log(alls)))
            else:
                s = sfn(x)
                
            for d in xrange(D):
                sup = sp.linspace(self.lb[d],self.ub[d],ns)
                X = sp.vstack([xmin]*ns)
                for j in xrange(ns):
                    X[j,d] = sup[j]
                [m,v] = self.G.infer_diag_post(X,[[sp.NaN]]*ns)
                
                sq = sp.sqrt(v)
                a[d].fill_between(sup,(m-2*sq).flatten(),(m+2.*sq).flatten(),facecolor = 'lightblue',edgecolor='lightblue')
                a[d].plot(sup,m.flatten())
                ps = sp.empty(ns)
                aps = sp.empty(ns)
                for j in xrange(ns):
                    ps[j] = self.query_pes(X[j,:],[s],[[sp.NaN]])
                    if d==0:
                        aps[j] = ps[j]/cfn(sp.array([X[j,:].flatten()]),s)
                
                a[d].twinx().plot(sup,ps,'r')
                if d==0:
                    a[d].twinx().plot(sup,aps,'g')
            
            for d in xrange(1,D):
                sup = sp.linspace(self.lb[d],self.ub[d],ns)
                X = sp.vstack([xmin]*ns)
                for j in xrange(ns):
                    X[j,d] = sup[j]
                    X[j,0] = 0.
                [m,v] = self.G.infer_diag_post(X,[[sp.NaN]]*ns)
                
                sq = sp.sqrt(v)
                a[d+D-1].fill_between(sup,(m-2*sq).flatten(),(m+2.*sq).flatten(),facecolor = 'lightblue',edgecolor='lightblue')
                a[d+D-1].plot(sup,m.flatten())
                ps = sp.empty(ns)
                aps = sp.empty(ns)
                for j in xrange(ns):
                    ps[j] = self.query_pes(X[j,:],[s],[[sp.NaN]])
                    
                
                a[d+D-1].twinx().plot(sup,ps,'r')
            
            
            a[2*D-1].hist(self.X[:,0].flatten(),50,facecolor='g')    
            print( xmin)
            
            f.savefig('../figcache/{0}.png'.format(time.time()))
            f.clf()
            plt.close(f)
            del(f)
            #plt.show()
        return [xmin,ymin,ierror]