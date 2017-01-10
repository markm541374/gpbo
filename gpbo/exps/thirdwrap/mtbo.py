from robo.fmin import mtbo
import numpy as np
import scipy as sp

import time
import sys
s=1e-6
tt0=time.time()
tc0=time.clock()

def f(x, **ev):
    # c = 1 - 0.5* ev['xa']
    c = 45 * sp.exp(-10. * ev['xa'])
    y = -sp.cos(x[0]) - sp.cos(x[1]) + 2. + 0.1 * s ** 2.
    b = ev['xa'] ** 2
    n = sp.random.normal() * sp.sqrt(s)
    if 'cheattrue' in ev.keys():
        if ev['cheattrue']:
            b = 0
            n = 0
    print 'f inputs x:{} ev:{} outputs y:{} (b:{} n:{}) c:{}'.format(x, ev, y + n + b, b, n, c)
    return y + b + n, c, dict()

def optmtbo(fn,lb,ub,salt,n,ninit=10,fname='results.csv'):
    D=len(ub)
    def objective_function(x, s):

        if s==1:
            y,c,aux = fn(x,**{'xa':0})
        elif s==0:
            y,c,aux = fn(x,**{'xa':salt})
        else:
            raise IndexError

        print "\nevaluation at {} {} returned {} {}\n".format(x,s,y,c)
        print "truetime {} clocktime {}".format(time.time()-tt0,time.clock()-tc0)
        sys.stdout.flush()
        sys.stderr.flush()
        return y, c

    res = mtbo(objective_function, lb, ub, n_tasks=2, n_init=ninit,num_iterations=n,burnin=100, chain_length=200)
    print res
    lf = open(fname, 'wb')
    lf.write(''.join(
        ['n, '] + ['x' + str(i) + ', ' for i in xrange(D)] + ['xa,']+ [
            'y, c, '] + ['rx' + str(i) + ', ' for i in xrange(D)] + [
            'truey at xrecc, taq, tev, trc, realtime']) + '\n')
    for i in xrange(n):
        st=''
        st+=str(i)+','
        st+='\n'
        lf.write(st)
    lf.close()

lb=sp.array([-1.,-1.])
ub=sp.array([ 1., 1.])
optmtbo(f,lb,ub,0.5,4,ninit=4)