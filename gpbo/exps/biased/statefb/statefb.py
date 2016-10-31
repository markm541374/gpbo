import scipy as sp
import time
import copy

def run(x0,n,T,con,trace=False):
    A = sp.array([[-0.5,0.],[0.,-1.]])
    B=sp.array([[1.],[1.]])
    C = sp.array([[con[0],con[1]]])
    t=0
    h=T/float(n)
    if trace:
        sys=sp.empty([3,n])
    x=sp.array([[x0[0]],[x0[1]]])
    r=0

    for i in xrange(n):
        #control signal
        u=-C.dot(x)
        #update the system
        dx = A.dot(x)+B.dot(u)

        x+=h*dx
        #state cost
        cs = x.T.dot(x)
        #control cost
        cc = 5.*u.dot(u.T)
        #cost update
        r+=(cs+cc)*h
        if trace:
            sys[:,i]=[x[0],x[1],u]
    if trace:
        return r[0,0],sys
    else:
        return r[0,0]

def f(x,**ev):
    l=-sp.log(2*10**-4)
    n=int((1e5)*sp.exp(-l*ev['xa']))
    t0=time.time()
    r = run([0.05,0.1],n,10,[2.*x[0],2.*x[1]],trace=False)
    t1=time.time()
    print 'running with C={} n={} (xa={})returned {} (10**{}) in {}'.format(x,n,ev['xa'],r,sp.log10(r),t1-t0)
    return sp.log10(r),t1-t0,dict()

def f_inplane(x,**ev):
    e = copy.deepcopy(ev)
    e['xa'] = 0
    y, c, aux = f(x, **e)
    return y, c, aux

def test():
    from matplotlib import pyplot as plt
    r,sys = run([0.05,0.1],10000,10,[0.05,0.3],trace=True)
    print r
    f,a = plt.subplots(2)
    a[0].plot(sys[0,:],'b')
    a[0].plot(sys[1,:],'r')
    a[1].plot(sys[2,:],'b')
    plt.show()

    return
def plot():
    import gpbo
    import os
    gpbo.core.ESutils.plot2dFtofile(f, os.path.join('dbout', 'trueobjective' + time.strftime(
        '%d_%m_%y_%H:%M:%S') + '.png'))
    return
if __name__=="__main__":
    test()
    #plot()


