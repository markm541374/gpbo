# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.
import os
import GPdc
import slice
import scipy as sp
from scipy import linalg as spl
from scipy import stats as sps
from matplotlib import pyplot as plt
from scipy.optimize import minimize as spomin
from scipy.stats import multivariate_normal as mnv
import time
import tqdm

SUPPORT_UNIFORM = 0
SUPPORT_SLICELCB = 1
SUPPORT_SLICEEI = 2
SUPPORT_SLICEPM = 3
SUPPORT_LAPAPR = 4
#drawing points between lb and ub using specified method
def draw_support(g, lb, ub, n, method, para=1.):
    
    #para is the std confidence bound
    if (type(g) is int):
        d=g
    else:
        d=g.D
    #print 'Draw support input GP with d={} lb {} ub {}'.format(d,lb,ub) 
    if method==SUPPORT_UNIFORM:
        print "Drawing support using uniform:"
        X=sp.random.uniform(size=[n,d])
        for i in xrange(d):
            X[:,i] *= ub[i]-lb[i]
            X[:,i] += lb[i]
    elif method==SUPPORT_LAPAPR:
        print "Drawing support using lapapr:"
        #start with 4 times as many points as needed
        #print 'a'
        para = 5*int(para)
        over = 4
        Xsto=sp.random.uniform(size=[over*para,d])
        for i in xrange(d):
            Xsto[:,i] *= ub[i]-lb[i]
            Xsto[:,i] += lb[i]
        #eval mean at the points
        #print 'b'
        fs = sp.empty(para*over)
        for i in xrange(para*over):
            fs[i] = g.infer_m_post(Xsto[i,:],[[sp.NaN]])[0,0]
        #print 'mean at random draws {}'.format(fs)
        Xst = sp.empty([2*para,d])
        #keep the lowest
        #print 'c'
        for i in xrange(para):
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
                print res
                if not res.status==2:
                    raise MJMError('failed in opt in support lapapr')
                else:
                    print "warn, lapapr opt did not fully vconverge "
            Xst[i+int(para),:] = res.x
            
        
        
        #find endpoints that are unique
        #print 'e'
        Xst
        unq = [Xst[0+para,:]]
        for i in xrange(para):
            tmp=[]
            for xm in unq:
                tmp.append(abs((xm-Xst[i+para,:])).max())
            
            if min(tmp)>0.0002:
                unq.append(Xst[i+para,:])
        #get the alligned gaussian approx of pmin
        #print 'f'
        #print unq
        cls = []
        for xm in unq:
            ls = []
            for i in xrange(d):
                
                vg = g.infer_diag_post(xm,[[i]])[1][0,0]
                gg = g.infer_diag_post(xm,[[i,i]])[0][0,0]
                ls.append(sp.sqrt(vg)/gg)
                
            cls.append(ls)
        #print 'g'
        X=mnv.rvs(size=n,mean=[0.]*d)
        if d==1:
            X.resize([n,1])
        neach = int(n/len(unq))
        for i in xrange(len(unq)):
            
            for j in xrange(d):
                X[i*neach:(i+1)*neach,j]*=min(2.,cls[i][j])
                X[i*neach:(i+1)*neach,j]+=unq[i][j]
            
        sp.clip(X,-1,1,out=X)
        if True:
        #if False:
            print "plotting draw_support...",
            if not os.path.exists(os.path.join('.', 'dbout')):
                os.mkdir('dbout')

            #print 'inits'
            #print Xst
            
            #print 'cls'
            #print cls
            
            np = para
            #print 'para{}'.format(para)
            #print Xst.shape
            n = 200
            x_ = sp.linspace(-1,1,n)
            y_ = sp.linspace(-1,1,n)
            z_ = sp.empty([n,n])
            s_ = sp.empty([n,n])
            for i in xrange(n):
                for j in xrange(n):
                    m_,v_ = g.infer_diag_post(sp.array([y_[j],x_[i]]),[[sp.NaN]])
                    z_[i,j] = m_[0,0]
                    s_[i,j] = sp.sqrt(v_[0,0])
            fig, ax = plt.subplots( nrows=2, ncols=1 ,figsize=(10,20))
            CS = ax[0].contour(x_,y_,z_,20)
            ax[0].clabel(CS, inline=1, fontsize=10)
            CS = ax[1].contour(x_,y_,s_,20)
            ax[0].axis([-1.,1.,-1.,1.])
            ax[1].clabel(CS, inline=1, fontsize=10)
            for i in xrange(np):
                ax[0].plot([Xst[i,0],Xst[i+np,0]],[Xst[i,1],Xst[i+np,1]],'b.-')
            for j in xrange(len(unq)):
                x = unq[j]
                ax[0].plot(x[0],x[1],'ro')
                xp = [x[0]+cls[j][0],x[0],x[0]-cls[j][0],x[0],x[0]+cls[j][0]]
                yp = [x[1],x[1]+cls[j][1],x[1],x[1]-cls[j][1],x[1]]
                ax[0].plot(xp,yp,'r-')
            fig.savefig(os.path.join('dbout','drawlapapr'+time.strftime('%d_%m_%y_%H:%M:%S')+'.png'))
            del(fig)
            print 'done'

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
        print "Drawing support using slice sample over LCB:"
        X = slice.slice_sample(f,0.5*(ub+lb),n,0.1*(ub-lb))
    
    elif method==SUPPORT_SLICEEI:
        def f(x):
            if all(x>lb) and all(x<ub):
                try:
                    ei=g.infer_EI(sp.array(x),[[sp.NaN]])[0,0]
                    #print ei
                    return sp.log(ei)
                except:
                    #ei=g.infer_EI(sp.array(x),[[sp.NaN]])[0,0]
                    g.printc()
                    raise
            else:
                return -1e99
        print "Drawing support using slice sample over EI:"
        X = slice.slice_sample(f,0.5*(ub+lb),n,0.1*(ub-lb))
    
    elif method==SUPPORT_SLICEPM:
        def f(x):
            if all(x>lb) and all(x<ub):
                [m,v] = g.infer_diag_post(sp.vstack([sp.array(x)]*d),[[i] for i in xrange(d)])
                p = 0.
                for i in xrange(d):
                    p+= -0.5*(m[0,i]**2)/v[0,i]
                ym = g.infer_m_post(sp.array(x),[[sp.NaN]])[0,0]
                if not sp.isfinite(p):
                    print [m,V,p]
                    #raise ValueError
                return -10*ym+0.01*p
            else:
                return -1e99
        if False:
            A = sp.empty([100,100])
            sup = sp.linspace(-0.999,0.999,100)
            for i in xrange(100):
                for j in xrange(100):
                    print sp.array([sup[i],sup[j]])
                    A[99-j,i] = f([sup[i],sup[j]])
                    print A[99-j,i]
            print A
            plt.figure()
            plt.imshow(A)
            plt.figure()
        print "Drawing support using slice sample over PM:"
        X = slice.slice_sample(f,0.5*(ub+lb),n,0.1*(ub-lb))
    else:
        raise RuntimeError("draw_support method invalid")
    return X


#return the min loc of draws on given support
def draw_min(g,support,n):
    print 'drawminin {}'.format(support)
    Z = g.draw_post(support, [[sp.NaN]]*support.shape[0],n)
    
    R = sp.empty([Z.shape[0],support.shape[1]])
    args = []
    for i in xrange(Z.shape[0]):
        a = sp.argmin(Z[i,:])
        args.append(a)
        R[i,:] = support[a,:]
    from itertools import groupby
    amins = [len(list(group)) for key, group in groupby(sorted(args))]
    print "In drawmin with {} support drew {} unique mins. Most freqent min chosen {}%".format(support.shape[0],len(amins),100.*max(amins)/float(n))
    

    if True:
    #if False:
        print 'plotting draw_min...',
        if not os.path.exists(os.path.join('.', 'dbout')):
            os.mkdir('dbout')
        from matplotlib import pyplot as plt
        import time
        #2d plot assuming [-1,1]^2 support
        n = 200
        x = sp.linspace(-1,1,n)
        y = sp.linspace(-1,1,n)
        z = sp.empty([n,n])
        s = sp.empty([n,n])
        for i in xrange(n):
            for j in xrange(n):
		m_,v_ = g.infer_diag_post(sp.array([0,y[j],x[i]]),[[sp.NaN]])
                z[i,j]=m_[0,0]
                s[i,j]=sp.sqrt(v_[0,0])
        fig, ax = plt.subplots( nrows=2, ncols=1 ,figsize=(10,20))
        ax[0].contour(x,y,z,20)
        CS = ax[1].contour(x,y,s,15)
        ax[1].clabel(CS, inline=1, fontsize=10)
        for i in xrange(Z.shape[0]):
            ax[0].plot(R[i,1],R[i,2],'ro')
        fig.savefig(os.path.join('dbout','drawmin'+time.strftime('%d_%m_%y_%H:%M:%S')+'.png'))
        del(fig)
        print 'done'
    return R

#fake gp class that 9looks like a d-1 gp becuase an extra vaue is added before callind
class gpfake():
    def __init__(self,g,axis,value):
        self.g=g
        self.D = g.D-1
        self.axis=axis
        self.value=value
        return
    
    def augx(self,x):
        ax = sp.hstack([x[:self.axis],sp.array([self.value]*(x.size/self.D)).T,x[self.axis:]])
        return ax
    
    def infer_m_post(self,x,d):
        return self.g.infer_m_post(self.augx(x),d)
    
    def infer_diag_post(self,x,d):
        return self.g.infer_diag_post(self.augx(x),d)
    
    def infer_EI(self,x,d):
        return self.g.infer_EI(self.augx(x),d)
    
    def infer_LCB_post(self,x,d,p):
        return self.g.infer_LCB_post(self.augx(x),d,p)
    
#ub and lb are still for the full space but the values in the chosen axis do not determine the outcome
def draw_support_inplane(g,lb,ub,n,method,axis,value,para=1.):
    print 'dsinplane axis:{} value:{}'.format(axis,value)
    
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
def drawhyp_plk(X,Y,S,D,ki,hm,hs,n,burn=80,subsam=5):
    def f(loghyp):
        ub = hm+1.8*hs
        lb = hm-1.8*hs
        if all(loghyp<ub) and all(loghyp>lb):
            r=GPdc.GP_LKonly(X,Y,S,D,GPdc.kernel(ki,X.shape[1],[10**i for i in loghyp])).plk(hm,hs)
            if sp.isnan(r):
                class MJMError(Exception):
                    pass
                print 'nan from GPLKonly with input'
                print [X,Y,S,D,ki,hm,hs,n,burn,subsam]
                raise MJMError('nan from GPLKonly with input')
        else:
            r=-1e99
        #print [loghyp, r]
        return r
    X = slice.slice_sample(f,hm,n,0.05*hs,burn=burn,subsam=subsam)
    return 10**X

#take a random draw of X points and draw Y from the specified kernel
def gen_dataset(nt,d,lb,ub,kindex,hyp,s=1e-9):
    X = draw_support(d, lb,ub,nt,SUPPORT_UNIFORM)
    D = [[sp.NaN]]*(nt)
    kf = GPdc.kernel(kindex,d,hyp)
    Kxx = GPdc.buildKsym_d(kf,X,D)
    Y = spl.cholesky(Kxx,lower=True)*sp.matrix(sps.norm.rvs(0,1.,nt)).T+sp.matrix(sps.norm.rvs(0,sp.sqrt(s),nt)).T
    S = sp.matrix([s]*nt).T
    return [X,Y,S,D]