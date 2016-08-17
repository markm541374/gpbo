# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.
import pyximport
from numpy import get_include
pyximport.install(setup_args={'include_dirs': get_include()})
import scipy as sp
from scipy import linalg as spl
import time
import GPdc
from matplotlib import pyplot as plt
import DIRECT
import PES
import ESutils
def cosines(x,s,d):
    x.resize([1,x.size])
    assert(d==[sp.NaN])
    
    u = 1.6*x[0,0]-0.5
    v = 1.6*x[0,1]-0.5
    f = 1.-(u**2 + v**2 -0.3*sp.cos(3*sp.pi*u)-0.3*sp.cos(3*sp.pi*v)+0.7)
    
    if s==0.:
        noise = 0.
    else:
        noise = sp.random.normal(scale=sp.sqrt(s))
    return [f +noise,1.]

def digitfn():
    from sklearn import datasets, svm
    digits = datasets.load_digits()
    n_samples = len(digits.images)
    data = digits.images.reshape((n_samples, -1))
    def digitojf(x,s,d,override=False):
        x.resize([1,x.size])
        g = 10**((x[0,0]-1.)*3)
        c = 6.*x[0,1]
        # Create a classifier: a support vector classifier
        classifier = svm.SVC(gamma=g, coef0=c,kernel='sigmoid')
        # We learn the digits on the first half of the digits
        classifier.fit(data[:n_samples / 2], digits.target[:n_samples / 2])
        # Now predict the value of the digit on the second half:
        r = classifier.score(data[n_samples / 2:],digits.target[n_samples / 2:])
        return [1-r,1.]
    return digitojf

def quad(x,s,d):
    assert(d==[sp.NaN])
    f = sum((x.flatten()-0.1)**2)
    if s==0.:
        noise = 0.
    else:
        noise = sp.random.normal(scale=sp.sqrt(s))
    return [f +noise,1.]

bananamin = sp.array([0.2,0.2])
def genbanana(ignores=-1.,cfn = lambda x:1.,):
    def banana(x,s,d, override=False):
        
        assert(d==[sp.NaN])
        x.resize([1,x.size])
        u = 5.*x[0,0]
        v = 5.*x[0,1]
        a=1.
        b=100.
        f = 1e-3*((a-u)**2 + b*(v-u**2)**2)
        if ignores>0:
            s=ignores
        if s==0.:
            noise = 0.
            
        else:
            
            noise = sp.random.normal(scale=sp.sqrt(s))
        try:
            return [f+noise,cfn(s)]
        except ZeroDivisionError:
            return [f+noise,sp.inf]
    return banana

def genbranin(ignores=-1.,cfn = lambda x:1.,):
    def branin(x,s,d, override=False):
        
        assert(d==[sp.NaN])
        x.resize([1,x.size])
        u = x[0,0]*7.5 + 2.5
        v = x[0,1]*7.5 + 2.5
        
        f = (-1.275*(u/sp.pi)**2+5*u/sp.pi+v-6)**2 +(10.-5./(4*sp.pi))*sp.cos(u) + 10.
        if ignores>0:
            s=ignores
        if s==0.:
            noise = 0.
            
        else:
            
            noise = sp.random.normal(scale=sp.sqrt(s))
        try:
            return [f+noise,cfn(s)]
        except ZeroDivisionError:
            return [f+noise,sp.inf]
    return branin

def gencamel(ignores=-1.,cfn = lambda x:1.,):
    def camel(x,s,d, override=False):
        
        assert(d==[sp.NaN])
        x.resize([1,x.size])
        u = x[0,0]*5
        v = x[0,1]*7.5
        
        f = 4*u**2+u*v-4*v**2-2.1*u**4+4*v**4+(u**6)/3.
        if ignores>0:
            s=ignores
        if s==0.:
            noise = 0.
            
        else:
            
            noise = sp.random.normal(scale=sp.sqrt(s))
        try:
            return [f+noise,cfn(s)]
        except ZeroDivisionError:
            return [f+noise,sp.inf]
    return camel


def gensquexpdraw(d,lb,ub,ignores=-1):
    nt=14
    [X,Y,S,D] = ESutils.gen_dataset(nt, d, lb, ub, GPdc.SQUEXP, sp.array([1.5] + [0.30] * d))
    G = GPdc.GPcore(X, Y, S, D, GPdc.kernel(GPdc.SQUEXP, d, sp.array([1.5] + [0.30] * d)))
    def obj(x,s,d,override=False):
        #print [x,s,d]
        if ignores>0:
            s=ignores
        if s==0. or override:
            noise = 0.
        else:
            noise = sp.random.normal(scale=sp.sqrt(s))
        print "EVAL WITH NOISE: "+str(noise) + "FROM S= "+str(s)
        return [G.infer_m(x,[d])[0,0]+noise,1.]
    def dirwrap(x,y):
        z = G.infer_m(x,[[sp.NaN]])[0,0]
        #z = obj(x,0.,[sp.NaN])
        return (z,0)
    [xmin,ymin,ierror] = DIRECT.solve(dirwrap,lb,ub,user_data=[], algmethod=1, maxf=89000, logfilename='/dev/null')
    
    return [obj,xmin,ymin]

def gensquexpIPdraw(d,lb,ub,sl,su,sfn,sls,cfn):
    #axis = 0 value = sl
    #d dimensional objective +1 for s
    nt=25
    #print sp.hstack([sp.array([[sl]]),lb])
    #print sp.hstack([sp.array([[su]]),ub])
    [X,Y,S,D] = ESutils.gen_dataset(nt, d + 1, sp.hstack([sp.array([[sl]]),lb]).flatten(), sp.hstack([sp.array([[su]]),ub]).flatten(), GPdc.SQUEXP, sp.array([1.5] + [sls] + [0.30] * d))
    G = GPdc.GPcore(X, Y, S, D, GPdc.kernel(GPdc.SQUEXP, d + 1, sp.array([1.5] + [sls] + [0.30] * d)))
    def obj(x,s,d,override=False):
        x = x.flatten()
        if sfn(x)==0. or override:
            noise = 0.
        else:
            noise = sp.random.normal(scale=sp.sqrt(sfn(x)))
        
        return [G.infer_m(x,[d])[0,0]+noise,cfn(x)]
    def dirwrap(x,y):
        z = obj(sp.array([[sl]+[i for i in x]]),sl,[sp.NaN],override=True)
        return (z,0)
    [xmin0,ymin0,ierror] = DIRECT.solve(dirwrap,lb,ub,user_data=[], algmethod=1, maxf=89000, logfilename='/dev/null')
    lb2 = xmin0-sp.ones(d)*1e-4
    ub2 = xmin0+sp.ones(d)*1e-4
    [xmin,ymin,ierror] = DIRECT.solve(dirwrap,lb2,ub2,user_data=[], algmethod=1, maxf=89000, logfilename='/dev/null')
    #print "RRRRR"+str([xmin0,xmin,ymin0,ymin,xmin0-xmin,ymin0-ymin])
    return [obj,xmin,ymin]




class opt(object):
    def __init__(self,objective,lb,ub,para=None,initstate=False):
        self.initstate=initstate
        self.ojf = objective
        self.lb = lb
        self.ub = ub
        self.d = lb.size
        try:
            self.d = para['d']
        except:
            pass
        self.X = sp.empty([0,self.d])
        self.Y = sp.empty([0,1])
        self.S = sp.empty([0,1])
        self.D = []
        
        self.R = sp.empty([0,self.d])
        self.C = []
        self.T = []
        self.Tr = []
        self.Ymin = []
        self.Xmin = sp.empty([0,self.d])
        self.Yreg = sp.empty([0,1])
        self.Rreg = sp.empty([0,1])
        self.sdefault = 1e-6
        self.init_search(para)
        return
    
    def setstate(self):
        self.X = self.initstate[0]
        self.Y = self.initstate[1]
        self.S = self.initstate[2]
        self.D = self.initstate[3]
        
        self.R = self.initstate[4]
        self.C = self.initstate[5]
        self.T = self.initstate[6]
        self.Tr = self.initstate[7]
        self.Ymin = self.initstate[8]
        self.Xmin = self.initstate[9]
        self.Yreg = self.initstate[10]
        self.Rreg = self.initstate[11]
        return
    
    def init_search(self,para):
        if self.initstate:
            self.setstate()
        else:
            self.method = "Default: RandomSearch"
        return
    
    def run_search(self):
        return self.run_search_random()
    
    def run_search_random(self):
        #print self.d
        
        #print self.ub
        #print self.lb
        
        xnext = sp.random.uniform(size=self.d)*(self.ub-self.lb)+self.lb
        snext = self.sdefault
        dnext = [sp.NaN]
        return [xnext,snext,dnext]
    
    def reccomend(self):
        return self.reccomend_random()
    
    def reccomend_random(self):
        i = sp.argmin(self.Y)
        return self.X[i,:]
    
    def query_ojf(self,x,s,d):
        return self.ojf(x,s,d)
    
    def step(self,random=False):
        t0=time.time()
        if random:
            [x,s,d] = self.run_search_random()
            print "random selected: "+str([x,s,d])
        else:
            [x,s,d] = self.run_search()
            print "search found: "+str([x,s,d])
        t1=time.time()
        if self.X.shape[0]==0:
            xr = (self.lb+self.ub)/2.
        elif random:
            xr = self.reccomend_random()
        else:
            xr = self.reccomend()
        print "reccomend point "+str(xr)
        t2=time.time()
        
        [y,c] = self.query_ojf(x,s,d)
        
        self.X = sp.vstack([self.X,x])
        self.Y = sp.vstack([self.Y,y])
        self.S = sp.vstack([self.S,s])
        self.D.append(d)
        
        self.R = sp.vstack([self.R,xr])
        self.Tr.append(t2-t1)
        self.C.append(c)
        self.T.append(t1-t0)
        
        self.Ymin.append(sp.amin(self.Y))
        
        self.Xmin = sp.vstack([self.Xmin,self.X[sp.argmin(self.Y),:]])
        self.Yreg = sp.vstack([self.Yreg,self.ojf(self.X[sp.argmin(self.Y),:],0.,[sp.NaN],override=True)[0]])
        self.Rreg = sp.vstack([self.Rreg,self.ojf(self.R[-1,:],0.,[sp.NaN],override=True)[0]])
        return

    def compX(self,xtrue):
        n = self.X.shape[0]
        R = sp.empty([2,n])
        for i in xrange(n):
            R[0,i] = spl.norm(self.X[i,:]-xtrue)
            R[1,i] = sp.amin(R[0,:i+1])
        return R
    def plot(self,truex,truey,ax,c):
        ax[0].plot(self.Ymin,c)
        ax[0].set_ylabel('Ymin')
        if not truex==None:
            n=self.X.shape[0]
            M = sp.empty(n)
            V = sp.empty(n)
            Z = sp.empty(n)
            for i in xrange(n):
                M[i] = spl.norm(self.X[i,:]-truex)
                V[i] = spl.norm(self.R[i,:]-truex)
                Z[i] = spl.norm(self.Xmin[i,:]-truex)
            ax[1].semilogy(M,c,label=str(type(self)))
            ax[1].set_ylabel('Xeval')
            #ax[1].legend(loc='upper center',ncol=2).draggable()
            ax[2].semilogy(V,c)
            ax[2].set_ylabel('Xrecc')
            ax[4].semilogy([sum(self.C[:i]) for i in xrange(n)],V,c)
            ax[4].set_ylabel('Xrecc/acccost')
            ax[2].semilogy(Z,c,linestyle='--')
            ax[2].set_ylabel('Xmin')
            ax[6].semilogx(sp.log10(self.Yreg-truey).flatten(),c,linestyle='--')
            ax[6].semilogx(sp.log10(self.Rreg-truey).flatten(),c)
            #print "recc:\n"+str((self.Yreg-truey).flatten())+"\nmmmmm\n"+str((self.Rreg-truey).flatten())
        
            ax[6].set_ylabel('Regret')
            ax[5].plot([sum(self.C[:i]) for i in xrange(n)],sp.log10(self.Yreg-truey).flatten(),c,linestyle='--')
            ax[5].plot([sum(self.C[:i]) for i in xrange(n)],sp.log10(self.Rreg-truey).flatten(),c)
            ax[5].set_ylabel('REgret/cost')
        
        ax[3].plot(self.C,c)
        ax[3].set_ylabel('cost')
        
         
        
        
        #print "recc:\n"+str((self.Yreg-truey).flatten())+"\nmmmmm\n"+str((self.Rreg-truey).flatten())
        
        return
    
class LCBMLE(opt):
    def init_search(self,para):
        self.kindex = para[0]
        self.mprior = para[1]
        self.sprior = para[2]
        self.volper = para[3]
        self.s = para[4]
        self.sdefault = para[4]
        ninit=para[5]
        if self.initstate:
            self.setstate()
        else:
            for i in xrange(ninit):
                self.step(random=True)
        return
    
    def run_search(self):
        
        MAP = GPdc.searchMAPhyp(self.X, self.Y, self.S, self.D, self.mprior, self.sprior, self.kindex)
        try:
            del(self.G)
        except:
            pass
        self.G = GPdc.GPcore(self.X, self.Y, self.S, self.D, GPdc.kernel(self.kindex, self.d, MAP))
        def directwrap(x,y):
            x.resize([1,self.d])
            
            a = self.G.infer_LCB(x,[[sp.NaN]],1.)[0,0]
            return (a,0)
        [xmin,ymin,ierror] = DIRECT.solve(directwrap,self.lb,self.ub,user_data=[], algmethod=1, volper=self.volper, logfilename='/dev/null')
        return [xmin,self.s,[sp.NaN]]
    def reccomend(self):
        def dirwrap(x,y):
            m  =self.G.infer_m(x,[[sp.NaN]])[0,0]
            return (m,0)
        [xmin,ymin,ierror] = DIRECT.solve(dirwrap,self.lb,self.ub,user_data=[], algmethod=1, volper=self.volper, logfilename='/dev/null')
        return xmin
    
class EIMLE(opt):
    def init_search(self,para):
        self.kindex = para[0]
        self.mprior = para[1]
        self.sprior = para[2]
        self.volper = para[3]
        self.s = para[4]
        self.sdefault = para[4]
        ninit=para[5]
        if self.initstate:
            self.setstate()
        else:
            for i in xrange(ninit):
                
                self.step(random=True)
                
        return
    
    def run_search(self):
        
        MAP = GPdc.searchMAPhyp(self.X, self.Y, self.S, self.D, self.mprior, self.sprior, self.kindex)
        try:
            del(self.G)
        except:
            pass
        self.G = GPdc.GPcore(self.X, self.Y, self.S, self.D, GPdc.kernel(self.kindex, self.d, MAP))
        
        def directwrap(x,y):
            x.resize([1,self.d])
            
            a = self.G.infer_lEI(x,[[sp.NaN]])
            #print [x,a]
            #print G.infer_diag_post(x,[[sp.NaN]])
            return (-a[0,0],0)
        [xmin,ymin,ierror] = DIRECT.solve(directwrap,self.lb,self.ub,user_data=[], algmethod=1,  volper = self.volper, logfilename='/dev/null')
        
        return [xmin,self.s,[sp.NaN]]
    
    def reccomend(self):
        def dirwrap(x,y):
            m  =self.G.infer_m(x,[[sp.NaN]])[0,0]
            return (m,0)
        [xmin,ymin,ierror] = DIRECT.solve(dirwrap,self.lb,self.ub,user_data=[], algmethod=1, volper= self.volper, logfilename='/dev/null')
        return xmin
    
class PESFS(opt):
    def init_search(self,para):
        self.para=para
        self.sdefault = para['s']
        if self.initstate:
            self.setstate()
        else:
            for i in xrange(para['ninit']):
                self.step(random=True)
        return
    
    def run_search(self):
        print "begin PESFS:"
        try:
            del(self.pesobj)
        except:
            pass
        self.pesobj = PES.PES(self.X,self.Y,self.S,self.D,self.lb.flatten(),self.ub.flatten(),self.para['kindex'],self.para['mprior'],self.para['sprior'],DH_SAMPLES=self.para['DH_SAMPLES'], DM_SAMPLES=self.para['DM_SAMPLES'], DM_SUPPORT=self.para['DM_SUPPORT'],DM_SLICELCBPARA=self.para['DM_SLICELCBPARA'],mode=self.para['SUPPORT_MODE'])
        [xmin,ymin,ierror] = self.pesobj.search_pes(self.sdefault,volper=self.para['volper'])
        return [xmin,self.para['s'],[sp.NaN]]
    
    def reccomend(self):
        def dirwrap(x,y):
            m  =self.pesobj.G.infer_m(x,[[sp.NaN]])[0,0]
            return (m,0)
        [xmin,ymin,ierror] = DIRECT.solve(dirwrap,self.lb,self.ub,user_data=[], algmethod=1, volper=self.para['volper'], logfilename='/dev/null')
        return xmin
    
    def pcs(self):
        self.pesobj = PES.PES(self.X,self.Y,self.S,self.D,self.lb.flatten(),self.ub.flatten(),self.para['kindex'],self.para['mprior'],self.para['sprior'],DH_SAMPLES=self.para['DH_SAMPLES'], DM_SAMPLES=self.para['DM_SAMPLES'], DM_SUPPORT=self.para['DM_SUPPORT'],DM_SLICELCBPARA=self.para['DM_SLICELCBPARA'],mode=self.para['SUPPORT_MODE'])
        xmin = self.reccomend()
        plt.figure(1)
        plt.plot(xmin[0],xmin[1],'r.')
        print xmin
        plt.figure(2)
        plt.subplot(4,1,1)
        ns = 6000
        sup = sp.linspace(-1,1,ns)
        for i in xrange(2):
            X = sp.vstack([xmin for k in xrange(ns)])
            print X.shape
            for j in xrange(ns):
                X[j,i] = sup[j]
            [m,v] = self.pesobj.G.infer_diag_post(X,[[sp.NaN]]*ns)
            s = sp.sqrt(v)
            plt.subplot(4,1,2*i+1)
            plt.fill_between(sup,(m-2*s).flatten(),(m+2*s).flatten(), facecolor='lightblue',edgecolor='lightblue')
            plt.plot(sup,m.flatten())
            [m,v] = self.pesobj.G.infer_diag_post(X,[[i]]*ns)
            s = sp.sqrt(v)
            plt.subplot(4,1,2*i+2)
            plt.fill_between(sup,(m-2*s).flatten(),(m+2*s).flatten(), facecolor='lightblue',edgecolor='lightblue')
            plt.plot(sup,m.flatten(),'r')
            p = sp.exp(-0.5*(m**2)/v)
            
            
            
            plt.twinx().plot(sup,p.flatten(),'g')
        return
    
class PESIS(PESFS):
    def init_search(self,para):
        self.para=para
        self.sdefault = -1
        if self.initstate:
            self.setstate()
        else:
            for i in xrange(para['ninit']):
                self.step(random=True)
        return
    
    def run_search(self):
        print "begin PESIS:"
        try:
            #print "aaa"
            del(self.pesobj)
            #print "bbb"
        except:
            print "ccc"
        self.pesobj = PES.PES(self.X,self.Y,self.S,self.D,self.lb.flatten(),self.ub.flatten(),self.para['kindex'],self.para['mprior'],self.para['sprior'],DH_SAMPLES=self.para['DH_SAMPLES'], DM_SAMPLES=self.para['DM_SAMPLES'], DM_SUPPORT=self.para['DM_SUPPORT'],DM_SLICELCBPARA=self.para['DM_SLICELCBPARA'],mode=self.para['SUPPORT_MODE'],noS=True)
        [xmin,ymin,ierror] = self.pesobj.search_pes(-1,volper=self.para['volper'])
        return [xmin,0.,[sp.NaN]]


class PESVS(opt):
    def init_search(self,para):
        self.para=para
        self.sdefault = para['s']
        if self.initstate:
            self.setstate()
        else:
            for i in xrange(para['ninit']):
                self.step(random=True)
        return
    
    def run_search(self):
        print "begin PES:"
        try:
            del(self.pesobj)
        except:
            pass
        self.pesobj = PES.PES(self.X,self.Y,self.S,self.D,self.lb.flatten(),self.ub.flatten(),self.para['kindex'],self.para['mprior'],self.para['sprior'],DH_SAMPLES=self.para['DH_SAMPLES'], DM_SAMPLES=self.para['DM_SAMPLES'], DM_SUPPORT=self.para['DM_SUPPORT'],DM_SLICELCBPARA=self.para['DM_SLICELCBPARA'],mode=self.para['SUPPORT_MODE'])
        [Qmin,ymin,ierror] = self.pesobj.search_acq(self.para['cfn'],self.para['logsl'],self.para['logsu'],volper=self.para['volper'])
        return [Qmin[:-1],10**Qmin[-1],[sp.NaN]]
    
    def reccomend(self):
        def dirwrap(x,y):
            m  =self.pesobj.G.infer_m(x,[[sp.NaN]])[0,0]
            return (m,0)
        [xmin,ymin,ierror] = DIRECT.solve(dirwrap,self.lb,self.ub,user_data=[], algmethod=1, volper=self.para['volper'], logfilename='/dev/null')
        return xmin
    
    def query_ojf(self,x,s,d):
        [y,c0] = self.ojf(x,s,d)
        c = self.para['cfn'](x,s)
        return [y,c]

#I'm defining that the augmented dimension must be rescaled to [0,1] with 0 as the true objective
class PESIP(opt):
    def init_search(self,para):
        self.para=para
        self.sdefault = para['s']
        self.lb = sp.hstack([sp.array([[para['sl']]]),self.lb])
        self.ub = sp.hstack([sp.array([[para['su']]]),self.ub])
        print self.lb
        print self.ub
        if self.initstate:
            self.setstate()
        else:
            for i in xrange(para['ninit']):
                self.step(random=True)
        return
    
    def run_search(self):
        print "begin PES:"
        try:
            del(self.pesobj)
        except:
            pass
        self.pesobj = PES.PES_inplane(self.X,self.Y,self.S,self.D,self.lb.flatten(),self.ub.flatten(),self.para['kindex'],self.para['mprior'],self.para['sprior'],self.para['axis'],self.para['value'],DH_SAMPLES=self.para['DH_SAMPLES'], DM_SAMPLES=self.para['DM_SAMPLES'], DM_SUPPORT=self.para['DM_SUPPORT'],DM_SLICELCBPARA=self.para['DM_SLICELCBPARA'],mode=self.para['SUPPORT_MODE'])
        self.train_costest()
        [Qmin,ymin,ierror] = self.pesobj.search_acq(self.costest,self.para['sfn'],volper=self.para['volper'])
        return [Qmin,self.para['sfn'](Qmin),[sp.NaN]]
    
    def reccomend(self):
        def dirwrap(x,y):
            m  =self.pesobj.G.infer_m(sp.hstack([self.para['sl'],x]),[[sp.NaN]])[0,0]
            return (m,0)
        
        [xmin,ymin,ierror] = DIRECT.solve(dirwrap,self.lb[:,1:],self.ub[:,1:],user_data=[], algmethod=1, volper=self.para['volper'], logfilename='/dev/null')
        return sp.hstack([sp.array(self.para['sl']),xmin])
    def reccomend_random(self):
        i = sp.argmin(self.Y)
        
        return sp.hstack([sp.array([self.para['sl']]),self.X[i,1:]])
    
    def query_ojf(self,x,s,d):
        [y,c0] = self.ojf(x,s,d)
        #c = self.para['cfn'](x,s)
        return [y,c0]
    def train_costest(self):
        print self.X[:,0]
        X = self.X[:,0].copy().reshape([len(self.C),1])
        C = sp.array(self.C)
        D = [[sp.NaN]]*len(self.C)
        S = sp.zeros([len(self.C),1])
        lbc = sp.array([-3.,-2.,-6.])
        ubc = sp.array([3.,2.,1.])
        MLEC =  GPdc.searchMLEhyp(X, sp.log(C), S, D, lbc, ubc, GPdc.SQUEXPCS, mx=10000)
        print "MLEC "+str(MLEC)
        self.ce = GPdc.GPcore(X.copy(), sp.log(C), S, D, GPdc.kernel(GPdc.SQUEXPCS, 1, sp.array(MLEC)))
        #self.ce = GPdc.GPcore(X,C,S,D,GPdc.kernel(GPdc.SQUEXPCS,1,sp.array([1.,0.3,1e-3])))
        #self.ce.printc()
        return
    def costest(self,x,s):
        m = sp.exp(self.ce.infer_m(sp.array(x[:,0]),[[sp.NaN]]))
        return m
    def plotcosts(self,a):
        #print self.X[:,0].flatten().shape
        #print sp.array(self.C).shape
        self.train_costest()
        a.plot(self.X[:,0].flatten(),sp.array(self.C).flatten(),'g.')
        sup = sp.linspace(0,1,100)
        mc,vc = self.ce.infer_diag(sup,[[sp.NaN]]*100)
        sc =sp.sqrt(vc)
        a.plot(sup,sp.exp(mc.flatten()),'b')
        a.fill_between(sup,sp.exp(mc.flatten()-2*sc).flatten(),sp.exp(mc.flatten()+2*sc).flatten(), facecolor='lightblue',edgecolor='lightblue')
        return
        
class PESIPS(PESIP):
    def init_search(self,para):
        self.para=para
        self.sdefault = -1
        self.lb = sp.hstack([sp.array([[para['sl']]]),self.lb])
        self.ub = sp.hstack([sp.array([[para['su']]]),self.ub])
        #print self.lb
        #print self.ub
        if self.initstate:
            self.setstate()
        else:
            for i in xrange(para['ninit']):
                self.step(random=True)
        return
    
    def run_search(self):
        print "begin PES:"
        try:
            del(self.pesobj)
        except:
            pass
        self.pesobj = PES.PES_inplane(self.X,self.Y,self.S,self.D,self.lb.flatten(),self.ub.flatten(),self.para['kindex'],self.para['mprior'],self.para['sprior'],self.para['axis'],self.para['value'],DH_SAMPLES=self.para['DH_SAMPLES'], DM_SAMPLES=self.para['DM_SAMPLES'], DM_SUPPORT=self.para['DM_SUPPORT'],DM_SLICELCBPARA=self.para['DM_SLICELCBPARA'],mode=self.para['SUPPORT_MODE'],noS=True)
        self.train_costest()
        [Qmin,ymin,ierror] = self.pesobj.search_acq(self.costest,self.para['sfn'],volper=self.para['volper'])
        return [Qmin,0.,[sp.NaN]]
    def plot(self,truex,truey):
        f,a = plt.subplots(3)
        #print self.X.shape
        a[0].plot(self.X[:,0].flatten(),'b')
        a[0].set_ylabel("augx")
        a[0].twinx().plot(self.C,'r')
        a[1].plot(sp.log10(self.Rreg.flatten()-truey))
        a[1].set_ylabel("regret")
        a[2].plot([sum(self.C[:j]) for j in xrange(len(self.C))],sp.log10(self.Rreg.flatten()-truey))
        a[2].set_ylabel("regret/c")
        return
def bounds(Xs,Ys,ns=100):
    #use a gp to infer mean and bounds on sets of x/y data that have diffent x
    #f,a = plt.subplots(2)
    #for i in xrange(len(Ys)):
    #    a[0].plot(Xs[i],Ys[i])
    
    X = sp.hstack(Xs)
    np = X.size
    Y = sp.hstack(Ys)
    X.resize([np,1])
    Y.resize([np,1])
    #a[1].plot(X,Y,'r.')
    np = X.size
    S = sp.zeros(np)
    D = [[sp.NaN]]*np
    ki = GPdc.MAT52CS
    mprior = sp.array([1.,2.,1.])
    sprior = sp.array([2.,2.,2.])
    #MAPH = GPdc.searchMAPhyp(X,Y,S,D,mprior,sprior, ki,mx=500)
    MAPH = sp.array([0.5,5.,0.3])
    g = GPdc.GPcore(X, Y, S, D, GPdc.kernel(ki, 1, MAPH))
    sup = sp.linspace(min(X),max(X),ns)
    [m,V] = g.infer_diag_post(sup,[[sp.NaN]]*ns)
    std = sp.sqrt(V+MAPH[2])
    #plt.fill_between(sup.flatten(),(m-std).flatten(),(m+std).flatten(),facecolor='lightblue',edgecolor='lightblue',alpha=0.5)
    #a[1].plot(sup,m.flatten(),'b')
    return [sup,m,std]

from scipy.interpolate import interp1d
    
def mergelines(x,y):
    minx = max([min(i) for i in x])
    maxx = min([max(i) for i in x])
    fs = []
    for i in xrange(len(x)):
        #print [x[i].shape,y[i].shape]
        
        fs.append(interp1d(x[i],y[i]))
    X = [i for i in sorted(sp.hstack([sp.array(j) for j in x])) if i<=maxx and i>=minx]
    np = len(X)
    X=sp.array(X)
    Y=sp.empty(np)
    ub=sp.empty(np)
    lb=sp.empty(np)
    for i in xrange(np):
        q = [j(X[i]) for j in fs]
        Y[i] = sp.mean(q)
        v = sp.var(q)
        ub[i] = Y[i]+2.*sp.sqrt(v)
        lb[i] = Y[i]-2.*sp.sqrt(v)
    
    return X,Y,lb,ub