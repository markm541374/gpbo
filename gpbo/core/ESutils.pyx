# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.
#cython: profile=True
from __future__ import print_function
xrange=range
import os
from gpbo.core import GPdc
from gpbo.core import slice
import scipy as sp
from scipy import linalg as spl
from scipy import stats as sps
import gpbo
import logging
logger = logging.getLogger(__name__)
try:
    from matplotlib import pyplot as plt
    plots=True
except ImportError:
    plots=False
    plt=None

from scipy.optimize import minimize as spomin
from scipy.stats import multivariate_normal as mnv
import time

import tqdm
from libc.math cimport log10, log, isnan
import copy
SUPPORT_UNIFORM = 0
SUPPORT_SLICELCB = 1
SUPPORT_SLICEEI = 2
SUPPORT_SLICEPM = 3
SUPPORT_LAPAPR = 4
SUPPORT_LAPAPROT = 5
SUPPORT_VARREJ = 6
#drawing points between lb and ub using specified method
def draw_support(g, lb, ub, n, method, para=1.):
    n = int(n)

    #para is the std confidence bound
    if (type(g) is int):
        d=int(g)
    else:
        d=int(g.D)
    #print 'Draw support input GP with d={} lb {} ub {}'.format(d,lb,ub)
    cdef int i,j
    if method==SUPPORT_UNIFORM:
        print( "Drawing support using uniform:")
        X=sp.random.uniform(size=[n,d])
        for i in range(d):
            X[:,i] *= ub[i]-lb[i]
            X[:,i] += lb[i]
    elif method==SUPPORT_VARREJ:
        print( "Drawing support using varreject:")
        batch=500
        out=[]
        while len(out)<n:
            X=sp.random.uniform(size=[batch,d])
            for i in range(d):
                X[:,i] *= ub[i]-lb[i]
                X[:,i] += lb[i]
            pa = sp.random.uniform(size=batch)
            m,v = g.infer_diag_post(X,[[sp.NaN]]*batch)
            for j in range(batch):
                if para*pa[j]<v[0,j]:
                    out.append(X[j,:])
        X=sp.vstack(out[:n])
    elif method==SUPPORT_LAPAPR:

        print( "Drawing support using lapapr:")
        #start with 4 times as many points as needed
        #print 'a'
        para = int(para)
        over = 8
        Xsto=sp.random.uniform(size=[over*para,d])
        for i in range(d):
            Xsto[:,i] *= ub[i]-lb[i]
            Xsto[:,i] += lb[i]
        #eval mean at the points
        #print 'b'
        fs = sp.empty(para*over)
        for i in range(para*over):
            fs[i] = g.infer_m_post(Xsto[i,:],[[sp.NaN]])[0,0]
        #print 'mean at random draws {}'.format(fs)
        Xst = sp.empty([2*para,d])
        #keep the lowest
        #print 'c'
        for i in range(para):
            j = fs.argmin()
            Xst[i,:] = Xsto[j,:]
            fs[j]=1e99
        #minimize the posterior mean from each start
        #print 'd'
        def f(x):
            y= g.infer_m_post(sp.array(x),[[sp.NaN]])[0,0]
            bound=0
            r = max(abs(x))
            if r>1:
                #print 'offedge {}'.format(x)
                bound=(1e3*(r-1))**6
            return y+bound
        for i in range(para):
            res = spomin(f,Xst[i,:],method='Nelder-Mead',options={'xtol':0.0001,'maxfev':2000})
            if not res.success:
                class MJMError(Exception):
                    pass
                print( res)
                if not res.status==2:
                    raise MJMError('failed in opt in support lapapr')
                else:
                    print( "warn, lapapr opt did not fully vconverge ")
            Xst[i+int(para),:] = res.x
            
        
        
        #find endpoints that are unique
        #print 'e'

        unq = [Xst[0+para,:]]
        for i in range(para):
            tmp=[]
            for xm in unq:
                tmp.append(abs((xm-Xst[i+para,:])).max())
            
            if min(tmp)>0.0002:
                unq.append(Xst[i+para,:])


        cls = []
        for xm in unq:
            ls = []
            for i in range(d):
                
                vg = g.infer_diag_post(xm,[[i]])[1][0,0]
                gg = g.infer_diag_post(xm,[[i,i]])[0][0,0]
                ls.append(sp.sqrt(vg)/gg)
                
            cls.append(ls)
        #print 'g'
        X=mnv.rvs(size=n,mean=[0.]*d)
        if d==1:
            X.resize([n,1])
        neach = int(n/(len(unq)+1))
        for i in range(len(unq)):
            
            for j in range(d):
                X[i*neach:(i+1)*neach,j]*=min(2.,cls[i][j])
                X[i*neach:(i+1)*neach,j]+=unq[i][j]
        X[len(unq)*neach:,:]=sp.random.uniform(size=[n-len(unq)*neach,d])
        for i in range(d):
            X[len(unq)*neach:,i] *= ub[i]-lb[i]
            X[len(unq)*neach:,i] += lb[i]
        sp.clip(X,-1,1,out=X)
        from gpbo.core import debugoutput

        if debugoutput['drawlap'] and plots:
            print( "plotting draw_support...",)
            np = para
            #print 'para{}'.format(para)
            #print Xst.shape
            n = 200
            x_ = sp.linspace(-1,1,n)
            y_ = sp.linspace(-1,1,n)
            z_ = sp.empty([n,n])
            s_ = sp.empty([n,n])
            for i in range(n):
                for j in range(n):
                    m_,v_ = g.infer_diag_post(sp.array([y_[j],x_[i]]),[[sp.NaN]])
                    z_[i,j] = m_[0,0]
                    s_[i,j] = sp.sqrt(v_[0,0])
            fig, ax = plt.subplots( nrows=2, ncols=1 ,figsize=(10,20))
            CS = ax[0].contour(x_,y_,z_,20)
            ax[0].clabel(CS, inline=1, fontsize=10)
            CS = ax[1].contour(x_,y_,s_,20)
            ax[0].axis([-1.,1.,-1.,1.])
            ax[1].clabel(CS, inline=1, fontsize=10)
            for i in range(np):
                ax[0].plot([Xst[i,0],Xst[i+np,0]],[Xst[i,1],Xst[i+np,1]],'b.-')
            for j in range(len(unq)):
                x = unq[j]
                ax[0].plot(x[0],x[1],'ro')
                xp = [x[0]+cls[j][0],x[0],x[0]-cls[j][0],x[0],x[0]+cls[j][0]]
                yp = [x[1],x[1]+cls[j][1],x[1],x[1]-cls[j][1],x[1]]
                ax[0].plot(xp,yp,'r-')


            fig.savefig(os.path.join(debugoutput['path'],'drawlapapr'+time.strftime('%d_%m_%y_%H:%M:%S')+'.png'))
            fig.clf()
            plt.close(fig)
            del(fig)
            print( 'done')
    elif method==SUPPORT_LAPAPROT:
        print( "Drawing support using lapaprot:")
        #start with 4 times as many points as needed
        #print 'a'
        para = int(para)
        over = 8
        Xsto=sp.random.uniform(size=[over*para,d])
        for i in range(d):
            Xsto[:,i] *= ub[i]-lb[i]
            Xsto[:,i] += lb[i]
        #eval mean at the points
        #print 'b'
        fs = sp.empty(para*over)
        for i in range(para*over):
            fs[i] = g.infer_m_post(Xsto[i,:],[[sp.NaN]])[0,0]
        #print 'mean at random draws {}'.format(fs)
        Xst = sp.empty([2*para,d])
        #keep the lowest
        #print 'c'
        for i in range(para):
            j = fs.argmin()
            Xst[i,:] = Xsto[j,:]
            fs[j]=1e99
        #minimize the posterior mean from each start
        #print 'd'
        def f(x):
            y= g.infer_m_post(sp.array(x),[[sp.NaN]])[0,0]
            bound=0
            r = max(abs(x))
            if r>1:
                #print 'offedge {}'.format(x)
                bound=(1e3*(r-1))**6
            return y+bound
        for i in range(para):
            retries=5
            for j in xrange(retries):
                res = spomin(f,Xst[i,:]+sp.random.uniform(-0.01,0.01),method='Nelder-Mead',options={'xtol':0.00001,'maxfev':4000})
                #bounds=tuple([(lb[j],ub[j]) for j in range(d)])
                #res = spomin(f ,Xst[i,:]+sp.random.uniform(-0.01,0.01),method='L-BFGS-B',bounds=bounds,options={'gtol':0.00001,'maxfun':300})
                if not res.success:
                    logger.warning('failed to find local mean {}'.format(res))
                else:
                    break
            if not res.success:
                class MJMError(Exception):
                    pass
                print( res)
                if not res.status==2:
                    raise MJMError('failed in opt in support lapapr')
                else:
                    print( "warn, lapapr opt did not fully vconverge ")
            Xst[i+int(para),:] = res.x

        #find endpoints that are unique
        #print 'e'

        unq = [Xst[0+para,:]]
        for i in range(para):
            tmp=[]
            for xm in unq:
                tmp.append(abs((xm-Xst[i+para,:])).max())

            if min(tmp)>0.0002:
                unq.append(Xst[i+para,:])

        U = []
        E = []
        for xm in unq:
            #gradient inference
            G,cG = g.infer_full_post(sp.vstack([xm]*d),[[i] for i in xrange(d)])
            divs = [[]]*int((d*(d+1)/2))
            #hessian inference
            k=0
            for i in xrange(d):
                for j in xrange(i+1):
                    divs[k]=[i,j]
                    k+=1
            vecH,cvecH = g.infer_full_post(sp.vstack([xm]*int(d*(d+1)/2)),divs)
            #build hesian matrix
            H = sp.empty(shape=[d,d])
            k=0
            for i in xrange(d):
                for j in xrange(i+1):
                    H[i,j]=H[j,i]=vecH[0,k]
                    k+=1
            #svd on cov of grad
            Ucg,Ecg,Vcg = spl.svd(cG)
            #new rotated covariance
            C0 = spl.solve(H,Ucg)
            varP = C0.dot(sp.diag(Ecg).dot(C0.T))
            #svd of new cov
            Uvp,Evp,Vvp=spl.svd(varP)
            U.append(Uvp)
            E.append(Evp)
            #print '\nat {}\ngrad\n{} \nvargrad\n{} \nhess\n{} \nsvdvargrad\n {}\n{}\nnewcov\n{}\newsvd\n{}\n{}'.format(xm,G,cG,H,Ucg,Ecg,varP,Uvp,Evp)


        #print 'g'

        X = sp.empty(shape=[n,d])
        if d==1:
            X.resize([n,1])
        neach = int(n/(len(unq)+1))
        for i in range(len(unq)):
            X[i*neach:(i+1)*neach,:]= U[i].dot(sp.diag(sp.sqrt(E[i])).dot(sp.random.normal(size=[d,neach]))).T
            for j in range(d):
                X[i*neach:(i+1)*neach,j]+=unq[i][j]
        X[len(unq)*neach:,:]=sp.random.uniform(size=[n-len(unq)*neach,d])
        for i in range(d):
            X[len(unq)*neach:,i] *= ub[i]-lb[i]
            X[len(unq)*neach:,i] += lb[i]

        sp.clip(X,-1,1,out=X)

        from gpbo.core import debugoutput

        if debugoutput['drawlap'] and plots:
            print( "plotting draw_support...",)

            np = para
            #print 'para{}'.format(para)
            #print Xst.shape
            n = 200
            x_ = sp.linspace(-1,1,n)
            y_ = sp.linspace(-1,1,n)
            z_ = sp.empty([n,n])
            s_ = sp.empty([n,n])
            for i in range(n):
                for j in range(n):
                    m_,v_ = g.infer_diag_post(sp.array([y_[j],x_[i]]),[[sp.NaN]])
                    z_[i,j] = m_[0,0]
                    s_[i,j] = sp.sqrt(v_[0,0])
            fig, ax = plt.subplots( nrows=2, ncols=1 ,figsize=(10,20))
            CS = ax[0].contour(x_,y_,z_,20)
            ax[0].clabel(CS, inline=1, fontsize=10)
            CS = ax[1].contour(x_,y_,s_,20)
            ax[0].axis([-1.,1.,-1.,1.])
            ax[1].clabel(CS, inline=1, fontsize=10)
            for i in range(np):
                ax[0].plot([Xst[i,0],Xst[i+np,0]],[Xst[i,1],Xst[i+np,1]],'b.-')

            circ=sp.empty([2,100])
            for j in range(len(unq)):
                x = unq[j]
                ax[0].plot(x[0],x[1],'ro')

                for i in xrange(100):
                    theta = 2.*sp.pi*i/99.
                    circ[:,i]=U[j].dot(sp.array([sp.sin(theta)*sp.sqrt(E[j][0]),sp.cos(theta)*sp.sqrt(E[j][1])]))+unq[j].T
                ax[0].plot(circ[0,:],circ[1,:],'r')



                #ax[0].plot([x[0],x[0]+(svd[j][2][0,0])*0.1],[x[1],x[1]+(svd[j][2][1,0])*0.1],'g')
                #ax[0].plot([x[0],x[0]+(svd[j][2][0,1])*0.1],[x[1],x[1]+(svd[j][2][1,1])*0.1],'g')
            fig.savefig(os.path.join(debugoutput['path'],'drawlapaprot'+time.strftime('%d_%m_%y_%H:%M:%S')+'.png'))
            fig.clf()
            plt.close(fig)
            del(fig)
            print( 'done')


    elif method==SUPPORT_SLICELCB:
        def f(x):
            if all(x>lb) and all(x<ub):
                try:
                    return -g.infer_LCB_post(sp.array(x),[[sp.NaN]],para)[0,0]
                except:
                    g.infer_LCB_post(sp.array(x),[[sp.NaN]],para)[0,0]
                    g.printc()
                    raise
            else:
                return -1e99
        print( "Drawing support using slice sample over LCB:")
        X = slice.slice_sample(f,0.5*(ub+lb),n,0.1*(ub-lb))
    
    elif method==SUPPORT_SLICEEI:
        def f(x):
            if all(x>lb) and all(x<ub):
                try:
                    ei=g.infer_EI(sp.array(x),[[sp.NaN]])[0,0]
                    if isnan(ei):
                        return 1e-99
                    #print ei
                    return log(ei)
                except:
                    #ei=g.infer_EI(sp.array(x),[[sp.NaN]])[0,0]
                    g.printc()
                    raise
            else:
                return -1e99
        print( "Drawing support using slice sample over EI:")
        X = slice.slice_sample(f,0.5*(ub+lb),n,0.125*(ub-lb))
    
    elif method==SUPPORT_SLICEPM:
        def f(x):
            if all(x>lb) and all(x<ub):
                [m,v] = g.infer_diag_post(sp.vstack([sp.array(x)]*d),[[i] for i in range(d)])
                p = 0.
                for i in range(d):
                    p+= -0.5*(m[0,i]**2)/v[0,i]
                ym = g.infer_m_post(sp.array(x),[[sp.NaN]])[0,0]
                if not sp.isfinite(p):
                    print( [m,p])
                    #raise ValueError
                return -10*ym+0.01*p
            else:
                return -1e99
        if False and plots:
            A = sp.empty([100,100])
            sup = sp.linspace(-0.999,0.999,100)
            for i in range(100):
                for j in range(100):
                    print( sp.array([sup[i],sup[j]]))
                    A[99-j,i] = f([sup[i],sup[j]])
                    print( A[99-j,i])
            print(A)
            plt.figure()
            plt.imshow(A)
            plt.figure()
        print( "Drawing support using slice sample over PM:")
        X = slice.slice_sample(f,0.5*(ub+lb),n,0.1*(ub-lb))
    else:
        raise RuntimeError("draw_support method invalid")
    return X


#return the min loc of draws on given support
def draw_min(g,support,n):
    #print 'drawminin {}'.format(support)
    Z = g.draw_post(support, [[sp.NaN]]*support.shape[0],n)
    
    R = sp.empty([Z.shape[0],support.shape[1]])
    args = []
    for i in range(Z.shape[0]):
        a = sp.argmin(Z[i,:])
        args.append(a)
        R[i,:] = support[a,:]
    from itertools import groupby
    amins = [len(list(group)) for key, group in groupby(sorted(args))]
    print( "In drawmin with {} support drew {} unique mins. Most freqent min chosen {}%".format(support.shape[0],len(amins),100.*max(amins)/float(n)))


    return R

def draw_min_xypairgrad(g,support,n,x):
    m,d=support.shape
    support = sp.vstack([support,[x]*(d+1)])
    Z = g.draw_post(support, [[sp.NaN]]*(m+1)+[[i] for i in range(d)],n)

    R = sp.empty([n,d])
    Y = sp.empty([n,2+d])
    args = []
    for i in range(n):
        a = sp.argmin(Z[i,:(m+1)])
        args.append(a)
        R[i,:] = support[a,:]
        Y[i,0] = Z[i,a]
        Y[i,1:] = Z[i,m:]
    from itertools import groupby
    amins = [len(list(group)) for key, group in groupby(sorted(args))]
    print( "In drawmin with {} support drew {} unique mins. Most freqent min chosen {}%".format(m,len(amins),100.*max(amins)/float(n)))
    #print R,Y
    return R,Y,args

#fake gp class that 9looks like a d-1 gp becuase an extra vaue is added before callind
class gpfake():
    def __init__(self,g,axis,value):
        self.g=g
        self.D = g.D-1
        self.axis=axis
        self.value=value
        #print 'GPFAKE axis{} value{}'.format(axis,value)
        return
    
    def augx(self,x):
        if x.ndim==1:
            ax = sp.hstack([x[:self.axis],sp.array([self.value]*int(x.size/self.D)).T,x[self.axis:]])
        else:
            ax = sp.hstack([x[:self.axis].reshape(x.shape[0],self.axis),sp.array([[self.value]*int(x.size/self.D)]).T,x[self.axis:].reshape(x.shape[0],self.D-self.axis)])
        return ax
    
    def infer_m_post(self,x,d):
        dstar = [[i+1 if i>=self.axis else i for i in e] for e in d]
        return self.g.infer_m_post(self.augx(x),dstar)
    
    def infer_diag_post(self,x,d):
        dstar = [[i+1 if i>=self.axis else i for i in e] for e in d]
        return self.g.infer_diag_post(self.augx(x),d)

    def infer_full_post(self,x,d):
        dstar = [[i+1 if i>=self.axis else i for i in e] for e in d]
        return self.g.infer_full_post(self.augx(x),d)

    def infer_EI(self,x,d):
        dstar = [[i+1 if i>=self.axis else i for i in e] for e in d]
        return self.g.infer_EI(self.augx(x),d)
    
    def infer_LCB_post(self,x,d,p):
        dstar = [[i+1 if i>=self.axis else i for i in e] for e in d]
        return self.g.infer_LCB_post(self.augx(x),d,p)
    
#ub and lb are still for the full space but the values in the chosen axis do not determine the outcome
def draw_support_inplane(g,lb,ub,n,method,axis,value,para=1.):
    #print 'dsinplane axis:{} value:{}'.format(axis,value)
    
    if (type(g) is int):
        gf = g-1
    else:
        gf = gpfake(g,axis,value)
        
    
    lb_red = sp.hstack([lb[:axis],lb[axis+1:]])
    ub_red = sp.hstack([ub[:axis],ub[axis+1:]])
    X = draw_support(gf,lb_red,ub_red,n,method,para=para)
    return sp.hstack([X[:,:axis],sp.ones([n,1])*value,X[:,axis:]])
    

def plot_gp(g,axis,x,d):
    [m,v] = g.infer_diag(x,d)
    s = sp.sqrt(v)
    axis.fill_between(x,sp.array(m-2.*s).flatten(),sp.array(m+2.*s).flatten(),facecolor='lightblue',edgecolor='lightblue')
    axis.plot(x,m.flatten(),'b')
    return 0

#draw hyperparameters given data from posterior likelihood
def drawhyp_plk(X,Y,S,D,ki,hm,hs,n,burn=80,subsam=5,chains=1,prior='lognorm'):
    if prior=='lognorm':
        ub = hm+2.8*hs
        lb = hm-2.8*hs
        def f(loghyp):
            cdef double i,r
            if all(loghyp<ub) and all(loghyp>lb):
                r=GPdc.GP_LKonly(X, Y, S, D, GPdc.kernel(ki, X.shape[1], [10 ** i for i in loghyp])).plk(hm, hs,shape=prior)
                if isnan(r):
                    class MJMError(Exception):
                        pass
                    print( 'nan from GPLKonly with input')
                    print( [X,Y,S,D,ki,hm,hs,n,burn,subsam])
                    raise MJMError('nan from GPLKonly with input')
            else:
                r=-1e99
            #print [loghyp, r]
            return r

        starts = sp.vstack([hm]*chains)
        print('using {} slice chains'.format(chains))
        for i in range(chains):
            starts[i,:]+=2*(sp.random.uniform(size=len(ub))-0.5)*hs
        X = sp.vstack([slice.slice_sample(f,starts[j,:],n/chains+1,0.05*hs,burn=burn,subsam=subsam) for j in xrange(chains)])
        return 10**X[:n,:]

    elif prior=='gamma':
        global count
        count=0
        ub = hm*hs+5*sp.sqrt(hm*hs**2)
        lb = sp.zeros_like(ub)
        def f(hyp):
            global count
            count = count+1
            cdef double i,r
            if all(hyp<ub) and all(hyp>lb):
                r=GPdc.GP_LKonly(X, Y, S, D, GPdc.kernel(ki, X.shape[1], hyp)).plk(hm, hs,shape=prior)
                if isnan(r):
                    class MJMError(Exception):
                        pass
                    print( 'nan from GPLKonly with input')
                    print( [X,Y,S,D,ki,hm,hs,n,burn,subsam])
                    raise MJMError('nan from GPLKonly with input')
            else:
                r=-1e99
            #print [loghyp, r]
            return r

        starts = sp.vstack([hm*hs]*chains)
        print('using {} slice chains'.format(chains))
        for i in range(chains):
            starts[i,:]+=0.25*(sp.random.uniform(size=len(ub))-0.5)*(hm*hs)
        X = sp.vstack([slice.slice_sample(f,starts[j,:],n/chains+1,0.05*sp.sqrt(hm*hs**2),burn=burn,subsam=subsam) for j in xrange(chains)])

        print('used {} plk evals'.format(count))
        return X[:n,:]


#take a random draw of X points and draw Y from the specified kernel
def gen_dataset(nt,d,lb,ub,kindex,hyp,s=1e-9):
    X = draw_support(d, lb,ub,nt,SUPPORT_UNIFORM)
    D = [[sp.NaN]]*(nt)
    kf = GPdc.kernel(kindex, d, hyp)
    Kxx = GPdc.buildKsym_d(kf, X, D,s=s)
    Y = spl.cholesky(Kxx,lower=True)*sp.matrix(sps.norm.rvs(0,1.,nt)).T+sp.matrix(sps.norm.rvs(0,sp.sqrt(s),nt)).T
    S = sp.matrix([s]*nt).T
    return [X,Y,S,D]

def plot2dFtofile(f,fname,xmin=[False],atxa=0.):
    if not plots:
        print( 'XXXplots disabled')
        return
    cdef int i,j
    n = 50
    x_ = sp.linspace(-1, 1, n)
    y_ = sp.linspace(-1, 1, n)
    z_ = sp.empty([n, n])
    s_ = sp.empty([n, n])
    ev={'s': 0, 'xa': atxa, 'd': [sp.NaN],'cheattrue':True}
    for i in xrange(n):
        for j in xrange(n):
            m_ = f(sp.array([y_[j], x_[i]]), **ev)
            z_[i, j] = m_[0]
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))
    CS = ax.contour(x_, y_, z_, 50)
    ax.clabel(CS, inline=1, fontsize=8)
    if any(xmin):
        ax.plot(xmin[0], xmin[1], 'ro')
    fig.savefig(fname)
    fig.clf()
    plt.close(fig)
    del (fig)
    return

def accumulate(x):
    y = copy.deepcopy(x)
    cdef int i
    for i in xrange(len(x)-1):
        y[i+1]+=y[i]
    return y


def quartsirregular(xdata,ydata,xtarget):
    n = len(xdata)
    inters=sp.empty(shape=[n,len(xtarget)])
    for i in xrange(n):
        inters[i,:]= sp.interp(xtarget,xdata[i].values.flatten(),ydata[i].values.flatten(),left=sp.NaN,right=sp.NaN)
    return sp.percentile(inters,25,axis=0),sp.percentile(inters,50,axis=0),sp.percentile(inters,75,axis=0)


def medianirregular(xdata,ydata,xtarget):
    return percentileirregular(xdata,ydata,xtarget,50)

def percentileirregular(xdata,ydata,xtarget,per):
    n = len(xdata)
    inters=sp.empty(shape=[n,len(xtarget)])
    for i in xrange(n):
        inters[i,:]= sp.interp(xtarget,xdata[i].values.flatten(),ydata[i].values.flatten(),left=sp.NaN,right=sp.NaN)
    return sp.percentile(inters,per,axis=0)