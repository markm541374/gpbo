import gpbo
import scipy as sp


#dimension
D=2
#random initialization
ninit=10
#const diagonal noise variance
s=1e-6

def f(x,**ev):

    y=-sp.cos(x[0])-sp.cos(x[1])+2
    c=1.
    n = sp.random.normal()*sp.sqrt(s)
    #this is used to check the true value at the imcumbent if available for plotting. Not counted in the optimization
    if 'cheattrue' in ev.keys():
        if ev['cheattrue']:
            n=0
    print 'f inputs x:{} ev:{} outputs y:{} (n:{}) c:{}'.format(x,ev,y+n,n,c)
    return y+n,c,dict()

nsteps=50
#pes with constant known noise
C=gpbo.core.config.pesfsdefault(f,D,nsteps,s,'results','pesfs.csv')
C.aqpara['nrandinit']=ninit
out = gpbo.search(C)
print out