# Use the log likelihood of hyperperameters as an objective function.


import ESutils

import GPdc
import scipy as sp
import numpy.random as npr
import time




class lkojf:
    def __init(self):
        return
    
    def makedata_default(self,n):
        d = 2
        ki = GPdc.SQUEXPCS
        hy = sp.array([1.2]+[0.25]*d+[1e-3])
        lb = sp.array([-1.]*d)
        ub = sp.array([1.]*d)
        self.makedata(n,d,ki,hy,lb,ub)
        return 
    
    def makedata(self,n,d,ki,hy,lb,ub):
        self.n=n
        self.d=d
        self.ki=ki
        self.hy=hy
        self.nhy = hy.size
        self.lb=lb
        self.ub=ub
        t0 = time.time()
        [X,Y,S,D] = ESutils.gen_dataset(n,d,lb,ub,ki,hy,s=hy[-1])
        t1 = time.time()
        print "Setuptime: "+str(t1-t0)
        S = sp.zeros(S.shape)
        self.X=X
        self.Y=Y
        self.S=S
        self.D=D
        return
        
    def llk(self,hy):
        t0 = time.clock()
        r = GPdc.GP_LKonly(self.X, self.Y, self.S, self.D, GPdc.kernel(self.ki, self.d, hy)).llk()
        t1 = time.clock()
        return [r,t1-t0]
    
    def llks(self,hy,sub):
        t0 = time.clock()
        Xs = sp.empty([sub,self.d])
        Ys = sp.empty([sub,1])
        Ss = sp.zeros([sub,1])
        Ds = []
        ix = npr.choice(range(self.n),size=sub, replace=False)
        for i,x in enumerate(ix):
            Xs[i,:] = self.X[x,:]
            Ys[i,:] = self.Y[x,:]
            Ds.append(self.D[x])
        t1 = time.clock()
        r = GPdc.GP_LKonly(Xs, Ys, Ss, Ds, GPdc.kernel(self.ki, self.d, hy)).llk()
        t2=time.clock()
        #print [t1-t0,t2-t1]
        return [r,t2-t0]