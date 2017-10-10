# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.
from __future__ import print_function
xrange=range
import scipy as sp
from scipy import linalg as spl
import logging
from gpbo.core import GPdc
import DIRECT
from scipy.optimize import minimize as spm
import gpbo
logger = logging.getLogger(__name__)


def trivialojf(x,**ev):
    return sum([(xi-0.1*i)**2 for i,xi in enumerate(x)]),1.,dict()

trivialymin = 0.

braninymin = 0.39788735772973816
def shiftbraninojf(x,**ev):
    if 'd' in ev.keys():
        assert(ev['d']==[sp.NaN])
    u = x[0]*7.5 + 2.5
    v = x[1]*7.5 + 2.5
        
    f = (-1.275*(u/sp.pi)**2+5*u/sp.pi+v-6)**2 +(10.-5./(4*sp.pi))*sp.cos(u) + 10.
    
    return f-braninymin,1.,dict()


hart3min = -3.86278
def shifthart3(x, **ev):
    #hartmann4 with a linear offset agains quadratic cost

    z = [0.5*xi +0.5 for xi in x]
    al = sp.array([1.,1.2,3.,3.2]).T
    A = sp.array([[3.,   10., 30.],
                  [0.1,  10., 35.],
                  [3.,   10., 30.],
                  [0.1,  10., 35.]])
    P = 0.0001*sp.array([[3689, 1170, 2673],
                         [4699, 4387, 7470],
                         [1091, 8732, 5547],
                         [381, 5743, 8828]])

    outer = 0
    for ii in range(4):
        inner = 0
        for jj in range(3):
            xj = z[jj]
            Aij = A[ii, jj]
            Pij = P[ii, jj]
            inner += Aij * (xj - Pij)**2



        new = al[ii] * sp.exp(-inner)
        outer = outer + new
    f = -outer

    print( 'f inputs x:{} ev:{} outputs y:{}'.format(z, ev, f))
    return f-hart3min, 1., dict()

def shifthart3D(x,**ev):
    f,c,d = shifthart3(x,**ev)
    h=1e-6
    F = [f]
    for i in range(3):
        zp = sp.copy(x)
        zp[i]+=h
        zn = sp.copy(x)
        zn[i]-=h
        print(zn,zp)
        fp,_,__=shifthart3(zp,**ev)
        fn,_,__=shifthart3(zn,**ev)
        F.append((fp-fn)/(2.*h))
    return F,c,d

hart6min = -3.32237
hart6xmin = [0.20169,0.150011,0.476874,0.275332,0.311652,0.6573]
def shifthart6(x, **ev):
    #hartmann4 with a linear offset agains quadratic cost
    z = [0.5*xi +0.5 for xi in x]
    al = sp.array([1.,1.2,3.,3.2]).T
    A = sp.array([[10., 3., 17., 3.5, 1.7, 8.],
                  [0.05, 10., 17., 0.1, 8, 14.],
                  [3., 3.5, 1.7, 10., 17., 8.],
                  [17., 8., 0.05, 10., 0.1, 14.]])
    P = 0.0001 * sp.array([[1312., 1696., 5569., 124., 8283., 5886.],
                           [2329., 4135., 8307., 3736., 1004., 9991.],
                           [2348., 1451., 3522., 2883., 3047., 6650.],
                           [4047., 8828., 8732., 5743., 1091., 381.]])

    outer = 0
    for ii in range(4):
        inner = 0
        for jj in range(6):
            xj = z[jj]
            Aij = A[ii, jj]
            Pij = P[ii, jj]
            inner += Aij * (xj - Pij)**2



        new = al[ii] * sp.exp(-inner)
        outer = outer + new
    f = -outer

    print( 'f inputs x:{} ev:{} outputs y:{}  '.format(z, ev, f))
    return f-hart6min, 1., dict()

def rosenojf(x,**ev):
    if 'd' in ev.keys():
        assert(ev['d']==[sp.NaN])
    if 's' in ev.keys():
        if ev['s']>0:
            noise = sp.random.normal(scale=sp.sqrt(ev['s']))
        else:
            noise=0
    else:
        noise=0.
    u = 5.*x[0]
    v = 5.*x[1]
    a=1.
    b=100.
    f = 1e-3*((a-u)**2 + b*(v-u**2)**2)
    return f+noise,1.,dict()

rosenxmin=[0.2,0.2]
rosenymin=0.

def genmat52ojf(d,lb,ub,A=1.,ls=0.3,fixs=-1,ki=GPdc.MAT52):
    from ESutils import gen_dataset
    if isinstance(ls,float):
        ls = [ls]*d
    #nt=sp.maximum(250,sp.minimum(int(d*100./sp.product(ls)),3000))
    nt=sp.maximum(150,sp.minimum(int(d*20./sp.product(ls)),2500))
    #nt=500

    s=1e-9

    for s in sp.logspace(-9,0,20):
        try:
            [X,Y,S,D] = gen_dataset(nt, d, lb, ub, ki, sp.array([A] + ls ),s=1e-9)
            break
        except:
            pass

   # from matplotlib import pyplot as plt
   # plt.figure()
   # plt.plot(sp.array(X),sp.array(Y),'b.')
   # plt.show(block=True)
    print('training GP')
    G = GPdc.GPcore(X, Y, S, D, GPdc.kernel(ki, d, sp.array([A] + ls )))

    def wrap(x):
        xq = sp.copy(x)
        xq.resize([1, d])
        a = G.infer_m_post(xq, [[sp.NaN]])
        return a[0, 0]

    dpara= {'user_data': [],
            'algmethod': 1,
            'maxf': 40000,
            'logfilename': '/dev/null'}
    lpara= {'ftol': 1e-20,
            'maxfun': 1200}
    xmin,ymin,ierror = gpbo.core.optutils.twopartopt(wrap,lb,ub,dpara,lpara)

    print('init {} {}'.format(xmin, ymin))
    for i in range(250):
        p = sp.random.normal(size=d)*1e-2
        res = spm( wrap,xmin+p,method='L-BFGS-B',bounds=tuple([(lb[j],ub[j]) for j in range(d)]),options={'gtol':1e-30})
       # print(res)
        #print(xmin,res.x,wrap(xmin),wrap(res.x)<wrap(xmin))
        if wrap(res.x) < wrap(xmin):
            xmin = res.x
        #    ymin = ymin+res.fun
            print('change: {} {}'.format(xmin,wrap(res.x)-ymin))
    ymin = wrap(xmin)
    def ojf(x,**ev):
        dx=ev['d']
        s=ev['s']
        if fixs<0:
            if ev['s']>0 and not 'cheattrue' in ev.keys():
                noise = sp.random.normal(scale=sp.sqrt(ev['s']))
            else:
                noise=0
        else:
            if not 'cheattrue' in ev.keys():
                noise = sp.random.normal(scale=sp.sqrt(fixs))
            else:
                noise=0.
        y= wrap(x)+noise#G.infer_m(sp.array(x),[dx])[0,0]+noise
        if not 'silent' in ev.keys():
            print('ojf at {} {} returned {} noise {}'.format([i for i in x],ev,y,noise))
        return y-ymin,1.,dict()
    logger.info('generated function xmin {} ymin {}(shifted to 0.) opt:{}'.format(xmin, ymin, ierror))
    return ojf,xmin,0.
    
def genbiasedmat52ojf(d,lb,ub,xls,sls):
    #s normalised to 0 exact, 1
    from ESutils import gen_dataset
    nt=30
    [X,Y,S,D] = gen_dataset(nt, d + 1, lb + [0], ub + [1], GPdc.MAT52, sp.array([1.5] + [xls] * d + [sls]))
    
    G = GPdc.GPcore(X, Y, S, D, GPdc.kernel(GPdc.MAT52, d + 1, sp.array([1.5] + [0.30] * d + [sls])))
    def ojf(x,**ev):
        #print "\nojfinput: {} : {}".format(x,ev)
        dx=ev['d']
        s=ev['s']
        if ev['s']>0:
            noise = sp.random.normal(scale=sp.sqrt(ev['s']))
        else:
            noise=0
        #print 'noise in ojf {}'.format(noise)
        xa=ev['xa']
        x = sp.array(x)
        xfull = sp.hstack([x,xa])
        return G.infer_m(xfull,[dx])[0,0]+noise,1.,dict()
    
    def dirwrap(x,y):
        z = G.infer_m(sp.hstack(sp.array(x)+[0.]),[[sp.NaN]])[0,0]
        return (z,0)
    [xmin,ymin,ierror] = DIRECT.solve(dirwrap,lb,ub,user_data=[], algmethod=1, maxf=40000, logfilename='/dev/null')
    #print [xmin, ymin]
    def spowrap(x):
        z = G.infer_m(sp.hstack(sp.array(x) + [0.]), [[sp.NaN]])[0, 0]
        return z
    y = spm(spowrap, xmin, method='l-bfgs-b',bounds=[(-1,1)]*d,options={'ftol':1e-15})
    xmin = y.x
    ymin = spowrap(y.x)
    #print [xmin,ymin]

    if sp.isnan(ymin):
        logger.warning('generator got nan optimizing objective. retrying...')
        return genbiasedmat52ojf(d,lb,ub,xls,sls)
    logger.info('generated function xmin {} ymin {} globopt:{} locopt:{}'.format(xmin, ymin, ierror, y.status))
    return ojf, xmin, ymin

def coswithbias():
    def ojf(x,**ev):
        xa=ev['xa']
        cost = 0.5+xa**2
        y = -sp.cos(x[0]*0.75/sp.pi)*sp.cos(x[1]*0.75/sp.pi)+0.05*xa**2
        return y,cost,dict()
    return ojf,sp.array([0,0]),-1.

def costfnwrap(ojfbase,cfn):
    def ojf(x,**ev):
        y,c0,ojfaux = ojfbase(x,**ev)
        c = cfn(x,**ev)
        return y,c,ojfaux
    return ojf




def gendecayingpositiveojf(d, lb, ub,sim):

    # s normalised to 0 exact, 1
    from ESutils import gen_dataset
    nt = 20
    cl=2.

    [X, Y, S, D] = gen_dataset(nt, d, lb, ub , GPdc.MAT52, sp.array([1] + [0.30] * d))
    G0 = GPdc.GPcore(X, Y, S, D, GPdc.kernel(GPdc.MAT52, d, sp.array([1] + [0.30] * d)))

    [X, Y, S, D] = gen_dataset(nt, d, lb, ub, GPdc.MAT52, sp.array([1] + [0.30] * d))
    G1 = GPdc.GPcore(X, Y, S, D, GPdc.kernel(GPdc.MAT52, d, sp.array([1] + [0.30] * d)))

    [X, Y, S, D] = gen_dataset(nt, d, lb, ub, GPdc.MAT52, sp.array([1] + [0.30] * d))
    G2 = GPdc.GPcore(X, Y, S, D, GPdc.kernel(GPdc.MAT52, d, sp.array([1] + [0.30] * d)))

    def p0(x):
        v= G0.infer_m(x, [[sp.NaN]])[0,0]
        y = v+1 if v>0 else sp.exp(v)
        return y

    def p1(x):
        v= G1.infer_m(x, [[sp.NaN]])[0, 0]
        y = v + 1 if v > 0 else sp.exp(v)
        return y

    def p2(x):
        v= G2.infer_m(x, [[sp.NaN]])[0, 0]
        y = v + 1 if v > 0 else sp.exp(v)
        return y

    def ojf(x, **ev):
        #print "ex: {} {}".format(x,ev['xa'])
        # print "\nojfinput: {} : {}".format(x,ev)
        dx = ev['d']
        s = ev['s']
        if ev['s'] > 0:
            noise = sp.random.normal(scale=sp.sqrt(ev['s']))
        else:
            noise = 0
        # print 'noise in ojf {}'.format(noise)
        xa = ev['xa']
        x = sp.array(x)


        y0=p0(x)
        y1=sim*p1(x)+y0
        l=p2(x)

        A = (y1-y0)/(sp.exp(l)-1)
        B=y0-A
        y = A*sp.exp(l*xa)+B

        c = sp.exp(-cl * xa)
        return y + noise, c, dict()

    def dirwrap(x, y):
        z = ojf(x,**{'d':[sp.NaN],'s':0,'xa':0})[0]
        return (z, 0)

    [xmin, ymin, ierror] = DIRECT.solve(dirwrap, lb, ub, user_data=[], algmethod=1, maxf=20000, logfilename='/dev/null')

    # print [xmin, ymin]
    def spowrap(x):
        z = ojf(x,**{'d':[sp.NaN],'s':0,'xa':0})[0]
        return z

    y = spm(spowrap, xmin, method='nelder-mead', options={'fatol': 1e-12})
    xmin = y.x
    ymin = spowrap(y.x)
    # print [xmin,ymin]
    logger.info('generated function xmin {} ymin {} yminisnan{} globopt:{} locopt:{}'.format(xmin, ymin, sp.isnan(ymin),ierror, y.status))
    if sp.isnan(ymin):
        logger.warning('generator got nan optimizing objective. retrying...')
        return gendecayingpositiveojf(d, lb, ub)
    return ojf, xmin, ymin
