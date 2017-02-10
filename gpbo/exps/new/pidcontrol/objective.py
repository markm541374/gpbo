from __future__ import print_function
import numpy as np
import scipy as sp
import time
import invp2 as invp
import sys

def f(x, **ev):
    #branin with a linear offset agains quadratic cost
    s = -9+10.*ev['xa']
    p = 10**x[0]
    i = 10**x[1]
    d = 10**x[2]
    t0=time.clock()
    y = np.log10(invp.runsim([100.,0.1,0.8,p,i,d],step=s,plot=False))
    t1=time.clock()
    c=10*(t1-t0)
    sys.stdout.flush()
    sys.stderr.flush()
    print( 'f inputs x:{} ev:{} outputs y:{}  c:{}'.format(x, ev, y, c))
    return y, c, dict()

if __name__=="__main__":
    #sys.exit(invp.runsim([100.,0.1,8.,9.3920548314116203, 0.010095206950107316, 6.1408677948580963],step=-10,plot=True))
    from matplotlib import pyplot as plt
    n = 30
    x_ = sp.linspace(-1.5, 2.2, n)
    y_ = sp.linspace(-1.5, 2.2, n)
    z_ = sp.empty([n, n])
    s_ = sp.empty([n, n])
    for i in range(n):
        print( i)
        for j in range(n):
            #m_ = invp.runsim([100.,0.1,8.,9.,(10**(x_[j])*10**(y_[i])),10**(y_[i])],step=-4,plot=False)
            #m_ = invp.runsim([100.,0.1,8.,10**(x_[j]),0.1,10**(y_[i])],step=-3,plot=False)
            #z_[i, j] = np.log10(m_)
            alpha=100.
            beta=0.1
            gamma=8.
            a=10**(x_[j])
            b = 0.1
            c = 10**(y_[i])

            z_[i,j] =np.log10(((gamma+gamma*a+alpha*c)*(alpha-0.7*9.81+beta*c+alpha*a+b*gamma)*(beta+beta*a+alpha*b)/(0.5*0.3+c*gamma)-(beta+beta*a+alpha*b)**2)/(gamma+gamma*a+alpha*c) - beta*b)
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(10, 20))
    CS = ax[0].contour(x_, y_, z_, 20)
    ax[0].clabel(CS, inline=1, fontsize=10)
    #CS = ax[1].contour(x_, y_, s_, 20)
    #ax[0].axis([-1., 1., -1., 1.])
    plt.show()
    #ax[1].clabel(CS, inline=1, fontsize=10)