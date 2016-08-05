# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.
import scipy as sp
import logging
import GPdc
import DIRECT
logger = logging.getLogger(__name__)


def trivialojf(x,**ev):
    return sum([(xi-0.1*i)**2 for i,xi in enumerate(x)]),1.,dict()

trivialymin = 0.

def braninojf(x,**ev):
    if 'd' in ev.keys():
        assert(ev['d']==[sp.NaN])
    if 's' in ev.keys():
        noise = sp.random.normal(scale=sp.sqrt(ev['s']))
    else:
        noise=0.
    u = x[0]*7.5 + 2.5
    v = x[1]*7.5 + 2.5
        
    f = (-1.275*(u/sp.pi)**2+5*u/sp.pi+v-6)**2 +(10.-5./(4*sp.pi))*sp.cos(u) + 10.
    
    return f+noise,1.,dict()

braninymin = 0.39788735772973816

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

def genmat52ojf(d,lb,ub):
    from ESutils import gen_dataset
    nt=14
    [X,Y,S,D] = gen_dataset(nt,d,lb,ub,GPdc.MAT52,sp.array([1.5]+[0.55]*d))
    G = GPdc.GPcore(X,Y,S,D,GPdc.kernel(GPdc.MAT52,d,sp.array([1.5]+[0.55]*d)))
    def ojf(x,**ev):
        dx=ev['d']
        s=ev['s']
        if ev['s']>0:
            noise = sp.random.normal(scale=sp.sqrt(ev['s']))
        else:
            noise=0
        return G.infer_m(sp.array(x),[dx])[0,0]+noise,1.,dict()
    def dirwrap(x,y):
        z = G.infer_m(x,[[sp.NaN]])[0,0]
        #z = obj(x,0.,[sp.NaN])
        return (z,0)
    [xmin,ymin,ierror] = DIRECT.solve(dirwrap,lb,ub,user_data=[], algmethod=1, maxf=89000, logfilename='/dev/null')
    logger.info('generated function xmin {} ymin {} {}'.format(xmin,ymin,ierror))
    
    return ojf,xmin,ymin
    
def genbiasedmat52ojf(d,lb,ub,sls):
    #s normalised to 0 exact, 1
    from ESutils import gen_dataset
    nt=20
    [X,Y,S,D] = gen_dataset(nt,d+1,lb+[0],ub+[1],GPdc.MAT52,sp.array([1.5]+[0.30]*d+[sls]))
    
    G = GPdc.GPcore(X,Y,S,D,GPdc.kernel(GPdc.MAT52,d+1,sp.array([1.5]+[0.30]*d+[sls])))
    def ojf(x,**ev):
        #print "\nojfinput: {} : {}".format(x,ev)
        dx=ev['d']
        s=ev['s']
        if ev['s']>0:
            noise = sp.random.normal(scale=sp.sqrt(ev['s']))
        else:
            noise=0
        xa=ev['xa']
        x = sp.array(x)
        xfull = sp.hstack([x,xa])
        return G.infer_m(xfull,[dx])[0,0]+noise,1.,dict()
    
    def dirwrap(x,y):
        z = G.infer_m(sp.hstack(sp.array(x)+[0.]),[[sp.NaN]])[0,0]
        return (z,0)
    [xmin,ymin,ierror] = DIRECT.solve(dirwrap,lb,ub,user_data=[], algmethod=1, maxf=89000, logfilename='/dev/null')
    logger.info('generated function xmin {} ymin {}'.format(xmin,ymin))
    return ojf, xmin, ymin
    
def costfnwrap(ojfbase,cfn):
    def ojf(x,**ev):
        y,c0,ojfaux = ojfbase(x,**ev)
        c = cfn(x,**ev)
        return y,c,ojfaux
    return ojf

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
def traincfn(x,c):
    n = x.size
    g = GPdc.GPcore(x,c,sp.array([1e-1]*n),[[sp.NaN]]*n,GPdc.kernel(GPdc.MAT52,1,[1.,0.2]))
    
    return cfnobj(g)