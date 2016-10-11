import gpbo
import scipy as sp


D=2
n=20
s=1e-6
def f(x,**ev):

    y=-sp.cos(x[0])-sp.cos(x[1])+2
    c=1.
    n = sp.random.normal()*1e-3
    if 'cheattrue' in ev.keys():
        if ev['cheattrue']:
            n=0
    print 'f inputs x:{} ev:{} outputs y:{} (n:{}) c:{}'.format(x,ev,y+n,n,c)
    return y+n,c,dict()


C=gpbo.core.config.pesfslearns(f,D,n,s,'results','pesfs.csv')

out = gpbo.search(C)
print out