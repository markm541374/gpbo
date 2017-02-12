from __future__ import print_function
import numpy as np
import scipy as sp
import time
import invp
import sys

def f(x, **ev):
    #branin with a linear offset agains quadratic cost
    s = -8+5.*ev['xa']
    rp = [np.log10(6.86),3.5]
    rd = [-1,3.]
    p = 10**(rp[0]+(rp[1]-rp[0])*(x[0]*0.5+0.5))
    d = 10**(rd[0]+(rd[1]-rd[0])*(x[2]*0.5+0.5))
    ri = [-2,np.log10((p*d-d*6.68)/(0.5*0.3))]
    i = 10**(ri[0]+(ri[1]-ri[0])*(x[1]*0.5+0.5))

    rpx = [1.,4.]
    px = 10**(rpx[0]+(rpx[1]-rpx[0])*(x[3]*0.5+0.5))
    rix = [0.,3.]
    ix = 10**(rix[0]+(rix[1]-rix[0])*(x[4]*0.5+0.5))
    rdx = [-1.,2.]
    dx = 10**(rdx[0]+(rdx[1]-rdx[0])*(x[5]*0.5+0.5))
    t0=time.clock()
    y = np.log10(invp.runsim([p,i,d,px,ix,dx],step=s,plot=False))
    t1=time.clock()
    c=100*(t1-t0)
    sys.stdout.flush()
    sys.stderr.flush()
    print( 'f inputs x:{} ev:{} outputs y:{}  c:{}'.format(x, ev, y, c))
    return y, c, dict()

if __name__=="__main__":
    #sys.exit(invp.runsim([100.,0.1,8.,9.3920548314116203, 0.010095206950107316, 6.1408677948580963],step=-10,plot=True))
    from matplotlib import pyplot as plt
    n = 30
    x_ = sp.linspace(-0.5, 3.2, n)
    y_ = sp.linspace(-2.5, 3.2, n)
    z_ = sp.empty([n, n])
    s_ = sp.empty([n, n])
    for i in range(n):
        print( i)
        for j in range(n):
            #m_ = invp.runsim([100.,0.1,8.,9.,(10**(x_[j])*10**(y_[i])),10**(y_[i])],step=-4,plot=False)
            m_ = invp.runsim([10**(x_[j]),10.5,10**(y_[i]),0.,0.,0.],step=-3,plot=False)
            z_[i, j] = np.log10(m_)
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(10, 20))
    CS = ax[0].contour(x_, y_, z_, 20)
    ax[0].clabel(CS, inline=1, fontsize=10)
    #CS = ax[1].contour(x_, y_, s_, 20)
    #ax[0].axis([-1., 1., -1., 1.])
    plt.show()
    #ax[1].clabel(CS, inline=1, fontsize=10)