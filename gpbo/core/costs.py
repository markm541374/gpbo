import scipy as sp
from scipy import linalg as spl
import logging
import GPdc
import os
import time
import DIRECT
from scipy.optimize import minimize as spm
import gpbo
logger = logging.getLogger(__name__)

def cf42(x,**ev):
    return 42.

def cfpower(A,p):
    def cf(x,**ev):
        s=ev['s']
        return A*s**(-p)
    return cf

def cfaexp(A,p):
    def cf(x,**ev):
        s=ev['xa']
        return A*sp.exp(-s*p)
    return cf
class cfnobj():
    def __init__(self,g):
        self.g=g
        return
    def __call__(self,x,**ev):
        xa =  ev['xa']
        return self.g.infer_m(sp.array([[xa]]),[[sp.NaN]])[0,0]

class logcfnobj():
    def __init__(self,g):
        self.g=g
        return
    def __call__(self,x,**ev):
        xa =  ev['xa']
        return sp.exp(self.g.infer_m(sp.array([[xa]]),[[sp.NaN]])[0,0])

class logcfnobjfull():
    def __init__(self,g):
        self.g=g
        return
    def __call__(self,x,**ev):
        x =  sp.hstack([x,ev['xa']])
        return sp.exp(self.g.infer_m(x,[[sp.NaN]])[0,0])


def traincfn1d(x,c):
    n = x.size
    g = GPdc.GPcore(x, c, sp.array([1e-1] * n), [[sp.NaN]] * n, GPdc.kernel(GPdc.MAT52, 1, [1., 0.2]))

    if gpbo.core.debugoutput and gpbo.core.debugoptions['cost1d']:
        print 'plotting cost1d...'
        from gpbo.core import debugpath
        import os
        if not os.path.exists(debugpath):
            os.mkdir(debugpath)
        import time
        from matplotlib import pyplot as plt
        f,a=plt.subplots(1)
        low = min(0,min(x))
        high = max(1,max(x))
        xaxis = sp.linspace(low,high,100)
        y,cy = g.infer_diag_post(xaxis,[[sp.NaN]]*100)

        a.plot(xaxis,y[0,:],'b')
        s = 2.*sp.sqrt(cy)
        u=sp.empty(100)
        l=sp.empty(100)
        for i in xrange(100):
            s = sp.sqrt(cy[0,i])
            u[i]=y[0,i]+2.*s
            l[i]=y[0,i]-2.*s
        a.fill_between(xaxis,l,u,facecolor='lightblue',edgecolor='lightblue',alpha=0.5)
        for i in xrange(n):
            a.plot(x[i],c[i],'r.')
        f.savefig(os.path.join(debugpath, 'cost1d' + time.strftime('%d_%m_%y_%H:%M:%S') + '.png'))
        del(f)
    return cfnobj(g)

def traincfn1dll(x,c):
    #cost modeled in a latent space c=exp(cl)
    n = x.size
    cl=sp.log(c)
    MAP = GPdc.searchMAPhyp(x, cl, sp.array([1e-3] * n), [[sp.NaN]] * n, sp.array([1.,0.,-1.]), sp.array([2.,2.,2.]), GPdc.MAT52CS)
    print 'MAPhyp in costfn {}'.format(MAP)
    g = GPdc.GPcore(x, cl, sp.array([1e-3] * n), [[sp.NaN]] * n, GPdc.kernel(GPdc.MAT52CS,1,MAP))

    if gpbo.core.debugoutput and gpbo.core.debugoptions['cost1d']:
        print 'plotting cost1d...'
        from gpbo.core import debugpath
        import os
        if not os.path.exists(debugpath):
            os.mkdir(debugpath)
        import time
        from matplotlib import pyplot as plt
        f,a=plt.subplots(2)
        low = min(0,min(x))
        high = max(1,max(x))
        xaxis = sp.linspace(low,high,100)
        y,cy = g.infer_diag_post(xaxis,[[sp.NaN]]*100)


        s = 2.*sp.sqrt(cy)
        u=sp.empty(100)
        l=sp.empty(100)
        for i in xrange(100):
            s = sp.sqrt(cy[0,i])
            u[i]=y[0,i]+2.*s
            l[i]=y[0,i]-2.*s
        a[0].plot(xaxis, y[0, :], 'b')
        a[0].fill_between(xaxis,l,u,facecolor='lightblue',edgecolor='lightblue',alpha=0.5)
        for i in xrange(n):
            a[0].plot(x[i],cl[i],'r.')
        a[0].set_ylabel('latent')

        a[1].plot(xaxis, sp.exp(y[0, :]), 'b')
        a[1].fill_between(xaxis, sp.exp(l), sp.exp(u), facecolor='lightblue', edgecolor='lightblue', alpha=0.5)
        for i in xrange(n):
            a[1].plot(x[i], c[i], 'r.')
        a[1].set_ylabel('out')

        f.savefig(os.path.join(debugpath, 'cost1d' + time.strftime('%d_%m_%y_%H:%M:%S') + '.png'))
        f.clf()
        plt.close(f)
        del(f)
    return logcfnobj(g)

def traincfnfull(x,c):
    #cost modeled in a latent space c=exp(cl)
    n,d = x.shape
    cl=sp.log(c)
    MAP = GPdc.searchMAPhyp(x, cl, sp.array([1e-6] * n), [[sp.NaN]] * n, sp.array([1.]+[-0.]*d ), sp.array([2.]*(d+1)), GPdc.MAT52)
    print 'MAPhyp in costfn {}'.format(MAP)
    g = GPdc.GPcore(x, cl, sp.array([1e-3] * n), [[sp.NaN]] * n, GPdc.kernel(GPdc.MAT52,1,MAP))
    if gpbo.core.debugoutput and gpbo.core.debugoptions['cost1d']:
        print 'plotting cost...'
        from gpbo.core import debugpath
        import os
        if not os.path.exists(debugpath):
            os.mkdir(debugpath)
        import time
        from matplotlib import pyplot as plt
        f,a=plt.subplots(d,2)

        n=60
        x=sp.linspace(-1,1,n)
        xa=sp.linspace(0,1,n)
        z=sp.empty([n,n])
        vz=sp.empty([n,n])
        for D in xrange(d-1):
            for i in xrange(n):
                for j in xrange(n):
                    q = sp.zeros([1, d])
                    q[0,0]=xa[j]
                    q[0,D+1]=x[i]
                    m,v = g.infer_diag(q,[[sp.NaN]])
                    z[i,j]=sp.exp(m[0,0])
                    vz[i,j]=sp.sqrt(v[0,0])
            try:
                CS=a[D,0].contour(xa,x,z,30)
                a[D, 0].clabel(CS, inline=1, fontsize=8)
            except ValueError:
                pass
            try:
                CS=a[D, 1].contour(xa, x, vz, 30)
                a[D, 1].clabel(CS, inline=1, fontsize=8)
            except ValueError:
                pass

        for i in xrange(n):
            for j in xrange(n):
                q = sp.zeros([1, d])
                q[0,0] = 0.5
                q[0, 1] = x[j]
                q[0, 2] = x[i]
                m, v = g.infer_diag(q, [[sp.NaN]])
                z[i, j] = sp.exp(m[0, 0])
                vz[i, j] = sp.sqrt(v[0, 0])
        try:
            CS = a[d-1, 0].contour(xa, x, z, 30)
            a[d-1, 0].clabel(CS, inline=1, fontsize=8)
        except ValueError:
            pass
        try:
            CS = a[d-1, 1].contour(xa, x, vz, 30)
            a[d-1, 1].clabel(CS, inline=1, fontsize=8)
        except ValueError:
            pass

        f.savefig(os.path.join(debugpath, 'cost1d' + time.strftime('%d_%m_%y_%H:%M:%S') + '.png'))
        f.clf()
        plt.close(f)
        del(f)
    return logcfnobjfull(g)


def predictive1d(x,c,t,ofs,C):

    cfbase=traincfn1dll(x,c)

    #make a polynomial fit for acquisition time
    T = t[ofs:]
    n = len(T)
    if n == 0:
        p=sp.array([10000 * sp.mean(t),0,0,0])
    elif n < 4:
        p = sp.array([sp.mean(T), 0, 0, 0])
    # TODO this bit properly
    else:
        X = sp.empty([n, 4])
        for i in xrange(n):
            X[i, 0] = 1
            X[i, 1] = float(i-n)
            X[i, 2] = float(i-n) ** 2
            X[i, 3] = float(i-n) ** 3
        p = spl.pinv(X).dot(T)
        if p[3] <= 0:
            p[:3] = spl.pinv(X[:, :3]).dot(T)
            p[3] = 0
            if p[2] <= 0:
                p[:2] = spl.pinv(X[:, :2]).dot(T)
                p[2] = 0
                if p[1] <= 0:
                    p = sp.array([sp.mean(T), 0, 0, 0])
    t_est = lambda x: p[0] + p[1] * x + p[2] * x ** 2 + p[3] * x ** 3

    def npred(cev):
        rts = sp.roots([p[3]/4.,p[2]/3.,p[1]/2.,p[0]+cev,-C])

        nmax = max([sp.real(i) for i in rts if sp.isreal(i)])
        cbar = p[0]+p[1]*nmax/2.+p[2]*(nmax**2)/3.+p[3]*(nmax**3)/4.
        return nmax,cbar

    def cfout(x,**ev):
        cev = cfbase(x,**ev)
        ecaq = npred(cev)[1]
        return cev+ecaq

    if gpbo.core.debugoutput and gpbo.core.debugoptions['taq']:
        from gpbo.core import debugpath
        from matplotlib import pyplot as plt
        print 'plotting taq...'
        f, a = plt.subplots(3)
        a[0].plot(t)


        xaxis = sp.linspace(ofs-0.5,len(t)+2)
        yaxis = map(t_est, xaxis-ofs-n)
        a[0].plot(xaxis,yaxis,'r')

        xaxis=sp.linspace(0,1,50)
        caxis=[cfbase(0.,**{'xa':i}) for i in xaxis]
        yaxis = map(npred,caxis)
        a[1].plot(xaxis,[i[0] for i in yaxis],'r')
        a[1].twinx().plot(xaxis, [i[1] for i in yaxis], 'b')


        xaxis=sp.linspace(0,1,100)
        cbase = map(lambda x:cfbase(None,**{'xa':x}),xaxis)
        cadj = map(lambda x:cfout(None,**{'xa':x}),xaxis)
        a[2].plot(xaxis,cbase,'b')
        a[2].plot(xaxis,cadj,'g')
        f.savefig(os.path.join(debugpath, 'taq' + time.strftime('%d_%m_%y_%H:%M:%S') + '.png'))
        f.clf()
        plt.close(f)
        del (f)
    return cfout