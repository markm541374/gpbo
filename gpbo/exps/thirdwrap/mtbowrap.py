from __future__ import print_function
xrange=range
from robo.fmin import mtbo
import numpy as np
import scipy as sp
import os
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
    print( 'f inputs x:{} ev:{} outputs y:{} (b:{} n:{}) c:{}'.format(x, ev, y + n + b, b, n, c))
    return y + b + n, c, dict()

def optmtbo(fn,lb,ub,salt,n,ninit=10,fname='results.csv',fpath='.',mod=False):
    D=len(ub)
    log=[]
    tinit=time.clock()
    def objective_function(x, s):
        t0=time.clock()
        if s==1:
            y,c,aux = fn(x,**{'xa':0})
        elif s==0:
            y,c,aux = fn(x,**{'xa':salt})
        else:
            raise IndexError
        t1=time.clock()
        print( "\nevaluation at {} {} returned {} {}\n".format(x,s,y,c))
        print( "truetime {} clocktime {}".format(time.time()-tt0,time.clock()-tc0))
        sys.stdout.flush()
        sys.stderr.flush()
        log.append({'x':x,'s':s,'y':y,'c':c,'t0':t0,'t1':t1})
        return y, c

    res = mtbo(objective_function, lb, ub, n_tasks=2, n_init=ninit,num_iterations=n,burnin=100, chain_length=200,mod=mod)
    print( 'results: {}'.format(os.path.join(fpath,fname)))
    lf = open(os.path.join(fpath,fname),'w')
    lf.write(''.join(
        ['n, '] + ['x' + str(i) + ', ' for i in xrange(D)] + ['xa,']+ [
            'y, c, '] + ['rx' + str(i) + ', ' for i in xrange(D)] + [
            'truey at xrecc, taq, tev, trc, realtime']) + '\n')
    for i in xrange(n):
        st=''
        st+=str(i)+','

        for j in xrange(D):
            st+=str(log[i]['x'][j])+','
        if log[i]['s']==0:
            st+=str(salt)+','
        else:
            st+=str(0.)+','
        st+=str(log[i]['y'])+','
        st+=str(log[i]['c'])+','
        for j in xrange(D):
            st+=str(res['trajectory'][i+1][j])+','
        st+=str(fn(res['trajectory'][i+1][:D],**{'xa':0,'cheattrue':True})[0])+','
        if i==0:
            st+=str(log[0]['t0']-tinit-res['recctimes'][0])+','
        else:
            st+=str(log[i]['t0']-log[i-1]['t1']-res['recctimes'][i])+','
        st+=str(log[i]['t1']-log[i]['t0'])+','
        st+=str(res['recctimes'][i])+','
        st+=time.strftime('%H:%M:%S  %d-%m-%y')
        st+='\n'
        lf.write(st)
    lf.close()

#lb=sp.array([-1.,-1.])
#ub=sp.array([ 1., 1.])
#optmtbo(f,lb,ub,0.5,14,ninit=10)