#example code to optimize using EI and maximum likelihood hyperparameters  with fixed noise and bias
import gpbo
import scipy as sp


#dimensionality
D=2
#noise variance
s=1e-6
#number of step to take
n=100

#define a simple 2d objective in x which also varies with respect to the environmental variable
def f(x,**ev):
    y=-sp.cos(x[0])-sp.cos(x[1])+2
    #fixed cost
    c=1.
    #noise
    n = sp.random.normal()*1e-3
    #we want to check the noiseless value when evaluating performance
    if 'cheattrue' in ev.keys():
        if ev['cheattrue']:
            n=0
    print('f inputs x:{} ev:{} outputs y:{} (n:{}) c:{}'.format(x,ev,y+n,n,c))
    return y+n,c,dict()

#arguments to generate default config are objective function, dimensionality, number of steps, noise variance, result directory and result filename
C=gpbo.core.config.eimledefault(f,D,n,s,'results','eimle.csv')
out = gpbo.search(C)
print(out)
