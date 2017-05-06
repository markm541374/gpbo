import gpbo
import scipy as sp


gpbo.core.debugoutput=True
gpbo.core.debugoptions={'datavis':True,'drawlap':False,'cost1d':False,'ctaq':False,'support':False,'adaptive':False,'logstate':False}

D=2
n=50
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


C=gpbo.core.config.pesfsdefault(f,D,n,s,'results','pesfs.csv')
#C.aqfn = gpbo.core.acquisitions.vmaxaq
out = gpbo.search(C)
print out