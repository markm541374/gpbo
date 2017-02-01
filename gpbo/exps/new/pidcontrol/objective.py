from __future__ import print_function
import numpy as np
import scipy as sp
import time
import invp
import sys

def f(x, **ev):
    #branin with a linear offset agains quadratic cost
    s = -9+10.*ev['xa']
    p = 10**(x[0]*1.5+1.5)
    i = (10**(1.5*x[1]-1.5))*p
    d = 10**(x[2])*p
    t0=time.clock()
    y = np.log10(invp.runsim([100.,0.1,0.8,p,i,d],step=s,plot=False))
    t1=time.clock()
    c=10*(t1-t0)
    sys.stdout.flush()
    sys.stderr.flush()
    print( 'f inputs x:{} ev:{} outputs y:{}  c:{}'.format(x, ev, y, c))
    return y, c, dict()

if __name__=="__main__":
    sys.exit(invp.runsim([100.,0.1,8.,8.3920548314116203, 0.010095206950107316, 2.1408677948580963],step=-10,plot=True))
