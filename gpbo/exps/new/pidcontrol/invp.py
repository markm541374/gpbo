import sys
import numpy as np
from scipy import integrate
import collections
from matplotlib import pyplot as plt
import time
5
System = collections.namedtuple('System','M m l g kap kai kad kxp kxi kxd')
def State(x=0.,dx=0.,a=0.,da=0.):
    return np.array([x,dx,a,da]).T

def dX(X,F,S):
    """
    get the acceleration of the system
    """
    F += (X[0]+10.)*0.25

    x  = X[0]
    dx = X[1]
    a  = X[2]
    da = X[3]
    a0 = S.M+S.m
    a1 = -S.m*S.l*np.cos(a)
    a2 = -np.cos(a)
    a3 = S.l
    b0 = F-S.m*S.l*(da**2)*np.sin(a)
    b1 = S.g*np.sin(a)
    
    d2x = (a3*b0-a1*b1)/(a0*a3-a1*a2)
    d2a = (-a2*b0+a0*b1)/(a0*a3-a1*a2)
    dx = dx
    da = da
    return np.array([dx,d2x,da,d2a]).T

class Control(object):
    def __init__(self,**para):
        for key in para:
            setattr(self,key,para[key])
        return
    def u(self,X,S):
        raise NotImplementedError

class NoCon(Control):
    def u(self,X,S):
        return -0.0

class PCon(Control):
    def u(self,X,S):
        e = X[2]-self.setA
        return -e*self.kP

class PIDCon(Control):
    def __init__(self,**para):
        super(PIDCon,self).__init__(**para)
        self.eI = 0. # integrated error
        self.eL = np.NaN # previous error

    def u(self,X,S):
        e = X[2]-self.setA
        
        self.eI = min(1.,max(-1.,self.eI+e*self.dt))
        eI = self.eI
        if np.isnan(self.eL):
            eD=0.
        else:
            eD = (e-self.eL)/self.dt
        self.eL=e
        return -(e*self.kP+eD*self.kD+eI*self.kI)

class CascadePIDCon(Control):
    def __init__(self,**para):
        super(CascadePIDCon,self).__init__(**para)
        self.eaI = 0. # outer controller
        self.eaL = np.NaN
        
        self.exI = 0. # outer controller
        self.exL = np.NaN
        self.Fl = 0. #prev output
        self.Al = 0. #prev angle target
        return

    def u(self,X,S):
        targeta = 0.
        ea = X[2]-targeta
         
        self.eaI = min(1.,max(-1.,self.eaI+ea*self.dt))
        eaI = self.eaI
        if np.isnan(self.eaL):
            eaD=0.
        else:
            eaD = (ea-self.eaL)/self.dt
        self.eaL=ea
        F = -(ea*self.kaP+eaD*self.kaD+eaI*self.kaI)+self.uinner(X,S)
        F = min(self.Fl+self.dt*400.,max(self.Fl-self.dt*400.,F))
        F =  min(20.,max(-20.,F))
        self.Fl=F
        return F
    def uinner(self,X,S):
        ex = -(X[0]-self.xtarget) # negative because x>0 requires a>0 to return
        
        self.exI = min(5.,max(-5.,self.exI+ex*self.dt))
        exI = self.exI
        if np.isnan(self.exL):
            exD=0.
        else:
            exD = (ex-self.exL)/self.dt
        self.exL=ex
        atarget =  -(ex*self.kxP+exD*self.kxD+exI*self.kxI)
        #atarget = min(self.Al+self.dt,max(self.Al-self.dt,atarget))
        #atarget =  min(np.pi/4.,max(-np.pi/4.,atarget))
        self.Al=atarget
        return atarget

def ES(x0,S,con,h,n):
    X = np.empty([x0.shape[0],n+1])
    X[:,0]=x0
    T = np.empty([n+1])
    T[0]=0.
    C = np.empty([n+1])
    C[0] = con.u(X[:,0],S)
    for i in range(n):
        X[:,i+1]=X[:,i]+h*dX(X[:,i],C[i],S)
        if X[2,i+1]>np.pi:
            X[2,i+1]-=2*np.pi
        if X[2,i+1]<-np.pi:
            X[2,i+1]+=2*np.pi
        T[i+1]=T[i]+h
        C[i+1] = con.u(X[:,i+1],S)
    return X,C,T

def MS(x0,S,con,h,n):
    X = np.empty([x0.shape[0],n+1])
    X[:,0]=x0
    T = np.empty([n+1])
    T[0]=0.
    C = np.empty([n+1])
    C[0] = con.u(X[:,0],S)
    for i in range(n):
        X[:,i+1]=X[:,i]+h*dX( X[:,i]+0.5*h*dX(X[:,i],C[i],S),C[i],S)
        if X[2,i+1]>np.pi:
            X[2,i+1]-=2*np.pi
        if X[2,i+1]<-np.pi:
            X[2,i+1]+=2*np.pi
        T[i+1]=T[i]+h
        C[i+1] = con.u(X[:,i+1],S)
    return X,C,T


def stepchange(controlpara,plot=False):
    s = System(M=0.5,m=0.2,l=0.3,g=9.8)
    x0 = State(x=1.,dx=0.0,a=-0.1,da=0.)
    Tmax=40.
    dt = 0.0001
    n = int(Tmax/dt)
#    c = NoCon(x=0.)
#    c = PIDCon(setA=-0.,kP=100.,kD=8.,kI=0.,dt=dt)
    kaP,kaI,kaD,kxP,kxI,kxD=controlpara
    c = CascadePIDCon(kaP=kaP, kaD=kaD, kaI=kaI, kxP=kxP, kxI=kxI, kxD=kxD,  xtarget=0., dt=dt)
    X,C,T = MS(x0,s,c,dt,n)

    for i in range(n+1):
        x=X[:,n-i]
        if abs(x[0])>0.01 or abs(x[1])>0.01 or abs(x[2])>0.01 or abs(x[3])>0.01:
            break
    settletime=T[n-i-1]
    abs_err = np.sum(np.abs(X[0,:n-i]))*dt
    #costfn = settletime + 0.1*abs_err*(max(0.,settletime-Tmax/2.))**2
    costfn=abs_err
    print 'Settle time: {} Abs Err: {} CostFn: {}'.format(settletime,abs_err,costfn)
    if plot:
        f,a = plt.subplots(5,sharex=True)
        a[0].plot(settletime,0,'ro')
        
        a[0].plot(T,X[0,:],'b')
        a[1].plot(T,X[1,:],'b')
        a[2].plot(T,X[2,:],'b')
        a[3].plot(T,X[3,:],'b')
        a[4].plot(T,C,'b')
    
        a[0].set_ylabel('x')
        a[1].set_ylabel('dx')
        a[2].set_ylabel('a')
        a[3].set_ylabel('da')
        a[4].set_ylabel('F')
        plt.show()
    return costfn
def spopt():
    def f(x):
        print x
        cp = [100.,0.1,8.,5.*np.exp(x[0]),1.*np.exp(x[1]),8.*np.exp(x[2])]
        print cp
        return stepchange(cp,plot=False)

    from scipy.optimize import minimize
    print minimize(f,[0.,0.,0.],method='Nelder-Mead',options={'maxfeval':100})

def dZ(t,X,S):
    """
    get the acceleration of the system
    """
    #F = -(X[0]+10.)*0.15-X[1]*0.1
    x  = X[0]
    dx = X[1]
    a  = (X[2]+np.pi)%(2*np.pi)-np.pi
    da = X[3]
    cx = X[4]
    ca = X[5]
    
    csa = -(S.kap*a+S.kai*ca+S.kad*da)
    csx = -(S.kxp*x+S.kxi*cx+S.kxd*dx)
    Ft=csa-csx 
    F = max(-20.,min(20.,Ft))-1.
    #F += -(X[0]+5.)*0.25 -X[1]*0.05
    a0 = S.M+S.m
    a1 = -S.m*S.l*np.cos(a)
    a2 = -np.cos(a)
    a3 = S.l
    b0 = F-S.m*S.l*(da**2)*np.sin(a)
    b1 = S.g*np.sin(a)
    d2x = (a3*b0-a1*b1)/(a0*a3-a1*a2)
    d2a = (-a2*b0+a0*b1)/(a0*a3-a1*a2)
    
    dca = 0. if (ca<-1. and a<0.) or (ca>1. and a>0.) else a
    dcx = 0. if (cx<-5. and x<0.) or (cx>5. and x>0.) else x

    dC = a**2+0.1*abs(da) + x**2 + 0.1*abs(dx)# + 0.*0.1*abs(F+1)#+x**2
    return [dx,d2x,da,d2a,dcx,dca,F,dC]



def runsim(controlpara,step=-6,plot=False):
#   K = [100.,0.1,8.,5., 1.,8. ]
#   K = [100.0, 0.1, 8.0,26.750950865287365, 2.396342520518441, 11.540708489014371] 
#   r = stepchange(K,plot=True)
    print 'in runsim controlpara={}'.format(controlpara)
    sys.stdout.flush()
    kap,kai,kad,kxp,kxi,kxd=controlpara
    s = System(M=0.5,m=0.2,l=0.3,g=9.8,kap=kap,kai=kai,kad=kad,kxp=kxp,kxi=kxi,kxd=kxd)
    x0 = [0.1,0.0,0.1,0.,0.,0.,0.,0.]
    O = integrate.ode(dZ).set_integrator('dopri5',atol=10**step,rtol=0.,verbosity=1,nsteps=10000)
    O.set_initial_value(x0,0.).set_f_params(s)
    Tmax=80.
    dt=0.025
    n=int(Tmax/dt)
    X = np.empty([8,n+1]) 
    X[:,0]=x0
    T = np.empty(n+1)
    T[0]=0.
    U = np.empty(n+1)
    U[0]=dZ(0.,x0,s)[6]
    i=1
    while O.successful() and i<n+1:
        O.integrate(O.t+dt)
        X[:,i]=O.y
        T[i]=O.t
        U[i]=dZ(O.t,O.y,s)[6]
        i+=1
    cost = X[7,-1]
    if not O.successful():
        cost=10**10
    print 'simresult={}'.format(cost)
    #print 'sqcost={} odetime={}'.format(cost,time.clock()-t0)
    if plot:
        f,a = plt.subplots(6,sharex=True)
        
        a[0].plot(T,X[0,:],'b')
        a[1].plot(T,X[1,:],'b')
        a[2].plot(T,(X[2,:]+np.pi)%(2*np.pi)-np.pi,'b')
        a[3].plot(T,X[3,:],'b')
        a[4].plot(T,U,'r') 
        a[5].plot(T,X[7,:],'g') 
        a[0].set_ylabel('x')
        a[1].set_ylabel('dx')
        a[2].set_ylabel('a')
        a[3].set_ylabel('da')
        plt.show()
    return cost

if __name__=="__main__":
    #sys.exit(spopt())
    #sys.exit(testfn())
    #sys.exit(runsim([100.,0.1,8.,25.,2.,11.],plot=True))
    sys.exit(runsim([1017.0344581161357, 3577.1169075731746, 65.698930381568232, 100.15805241836431, 6.5595255493927871, 0.10031635464239824],plot=True))
