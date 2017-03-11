from __future__ import print_function
import gpbo.core.objectives as objectives
import scipy as sp
D=2
lb = [-1,-1]
ub=[1.,1.]
def genf(xls,els,cls,cA):
    #xls : function lengthscale
    #els : environemntal variable lenghtscale
    #cls : costfunction decay constant
    #cA : costfunction value at full
    ojfw, xmin, ymin = objectives.genbiasedmat52ojf(D,lb,ub,xls,els)

    def f(x, **ev):
        #branin with a linear offset agains quadratic cost
        y,_,__ = ojfw(x,**ev)
        y-=ymin
        c = cA*sp.exp(-cls*ev['xa'])
        print( 'f inputs x:{} ev:{} outputs y:{}  c:{}'.format(x, ev, y, c))
        return y, c, dict()
    return f

truemin = 0.
