import gpbo
import scipy as sp

def cfn(x,s):
    return s**-0.5

def f(x,**ev):

    y = sp.sin(6.*x[0])+x[0]*2-0.5*x[0]-sp.cos(x[1]*5-0.2*x[0])+sp.exp(x[1]-x[0])
    n = sp.random.normal()*sp.sqrt(ev['s'])
    #this is used to check the true value at the imcumbent if available for plotting. Not counted in the optimization
    if 'cheattrue' in ev.keys():
        if ev['cheattrue']:
            n=0
    c=cfn(x,ev['s'])
    print 'f inputs x:{} ev:{} outputs y:{} (n:{}) c:{}'.format(x,ev,y+n,n,c)
    return y+n,c,dict()

D=2
n=30
lsl=-10
lsu=-3
#pes with variable fidelity but unbiased noise
C=gpbo.core.config.pesvsdefault(f,cfn,D,n,lsl,lsu,'results','pesvs.csv')
out = gpbo.search(C)
print out
