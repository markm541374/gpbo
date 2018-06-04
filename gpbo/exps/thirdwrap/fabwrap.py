from __future__ import print_function
xrange=range
from robo.fmin import fabolas, fabolas_mod
import numpy as np
import scipy as sp
import os
import time
import sys
s=1e-6
tt0=time.time()
tc0=time.clock()

def f(x, **ev):
    sn=ev['xa']
    y = 3.-np.cos(2.9*(1-0.2*sn)*(x[0]-sn*0.5))-np.cos(3.1*(x[1]+sn*0.2))+sn**2
    c=0.5+3.*(1.-sn)**2

    print( 'f inputs x:{} ev:{} outputs y:{}  c:{}'.format(x, ev, y ,c))
    return y , c, dict()

def optfabolas(fn,lb,ub,n,ninit,fname='results.csv',fpath='.',mod=False,switchestimator=False,switchkernel=False,timelimit=np.Inf):
    D=len(ub)
    log=[]
    tinit=time.clock()
    def objective_function(x, s):
        t0=time.clock()
        #adjust s to 0-1 where 0 is full
        #sn=1-np.log10(s)/4.+0.5
        sn = 1-(s-100.)/(1000000.-100.)
        y,c,aux = fn(x,**{'xa':sn})

        t1=time.clock()
        print( "\nevaluation at {} {} returned {} {}\n".format(x,s,y,c))
        print( "truetime {} clocktime {}".format(time.time()-tt0,time.clock()-tc0))
        sys.stdout.flush()
        sys.stderr.flush()
        log.append({'x':x,'s':sn,'y':y,'c':c,'t0':t0,'t1':t1})
        return y, c

    s_min = 100
    s_max = 1000000
    if mod:
        res = fabolas_mod(objective_function, lower=lb, upper=ub, s_min=s_min,s_max=s_max,n_init=ninit,num_iterations=n,switchestimator=switchestimator,switchkernel=switchkernel,timelimit=timelimit)
    else:
        res = fabolas(objective_function, lower=lb, upper=ub, s_min=s_min,s_max=s_max,n_init=ninit,num_iterations=n)
    print( res)
    print( 'results: {}'.format(os.path.join(fpath,fname)))
    lf = open(os.path.join(fpath,fname),'w')
    lf.write(''.join(
        ['n, '] + ['x' + str(i) + ', ' for i in xrange(D)] + ['xa,']+ [
            'y, c, '] + ['rx' + str(i) + ', ' for i in xrange(D)] + [
            'truey at xrecc, taq, tev, trc, realtime']) + '\n')
    print( n)
    inclast = -999
    ylast = -9999
    if not 'recctimes' in res.keys():
        res['recctimes']=[0.]*n
    for i in xrange(len(log)):
        st=''
        st+=str(i)+','

        for j in xrange(D):
            st+=str(log[i]['x'][j])+','
        st+=str(log[i]['s'])+','
        st+=str(log[i]['y'])+','
        st+=str(log[i]['c'])+','
        for j in xrange(D):
            st+=str(res['incumbents'][i+1][j])+','
        if inclast!=list(res['incumbents'][i+1][:D]):
            print('eval incumbent')
            y = fn(res['incumbents'][i+1][:D],**{'xa':0,'cheattrue':True})[0]
            ylast=y
            inclast = list(res['incumbents'][i+1][:D])
        else:
            print('no change to incumbent')
            y=ylast
        st+=str(y)+','
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

def optfabolas_mod(fn, lb, ub, n, ninit, fname='results.csv', fpath='.',switchkernel=False,switchestimator=False,timelimit=np.Inf):
    return optfabolas(fn, lb, ub, n, ninit, fname=fname, fpath=fpath,mod=True,switchestimator=switchestimator,switchkernel=switchkernel,timelimit=timelimit)

#lb=sp.array([-1.,-1.])
#ub=sp.array([ 1., 1.])
#optfabolas(f,lb,ub,50,ninit=20,mod=True)