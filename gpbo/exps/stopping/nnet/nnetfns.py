from subprocess import call
import numpy as np
import re
import time

def linmap(x,lb,ub):
    return lb+(ub-lb)*0.5*(x+1)

def logmap(x,lb,ub):
    return np.exp(linmap(x,np.log(lb),np.log(ub)))


def nnet_cancer(xraw,**ev):
    x=np.empty(4)
    # nnsize \in [1 100]
    x[0] = logmap(xraw[0],1,100)
    # mu \in [0.000001, 100]
    x[1] = logmap(xraw[1],0.000001,100)
    # mu_dec \in [0.01, 1]
    x[2] = logmap(xraw[2],0.01,1)
    # mu_inc \in [1.01, 20]
    x[3] = logmap(xraw[3],1.01,20)
    print('nnet cancer. GPinput {} NNinput {}'.format(xraw,x))
    t0=time.clock()
    print('Matlab call...')
    S = call(' matlab  -nodesktop -nojvm -nosplash -r \'nnet_cancer({},{},{},{})\''.format(*x),shell=True,stdout=open('tmp.txt','w'))
    t1=time.clock()
    print('...Matlab return')
    L = open('tmp.txt','r').read()
    R = 1-float(re.search('XXX(.*)XXX',L).group(1))
    T = t1-t0
    print('Result {} Time {}'.format(R,T))
    return R,T,dict()

def nnet_boston(xraw,**ev):
    x=np.empty(4)
    # nnsize \in [1 100]
    x[0] = linmap(xraw[0],1,100)
    # mu \in [0.000001, 100]
    x[1] = linmap(xraw[1],0.000001,100)
    # mu_dec \in [0.01, 1]
    x[2] = linmap(xraw[2],0.01,1)
    # mu_inc \in [1.01, 20]
    x[3] = linmap(xraw[3],1.01,20)
    print('nnet boston. GPinput {} NNinput {}'.format(xraw,x))
    t0=time.clock()
    print('Matlab call...')
    S = call(' matlab  -nodesktop -nojvm -nosplash -r \'nnet_boston({},{},{},{})\''.format(*x),shell=True,stdout=open('tmp.txt','w'))
    t1=time.clock()
    print('...Matlab return')
    L = open('tmp.txt','r').read()
    R = -float(re.search('XXX(.*)XXX',L).group(1))
    T = t1-t0
    print('Result {} Time {}'.format(R,T))
    return R,T,dict()

#print(nnet_boston([0.,0.,0.,0.]))