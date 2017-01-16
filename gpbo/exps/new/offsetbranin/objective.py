from __future__ import print_function
import scipy as sp



def f(x, **ev):
    #branin with a linear offset agains quadratic cost
    sn = ev['xa']
    u = x[0] * 7.5 + 2.5 + 0.2*sn
    v = x[1] * 7.5 + 2.5 - 0.15*sn

    f = (-1.275 * (u / sp.pi) ** 2 + 5 * u / sp.pi + v - 6) ** 2 + (10. - 5. / (4 * sp.pi)) * sp.cos(u) + 10. + 0.5*sn


    c = 1. + 24. * (1. - sn) ** 2

    print( 'f inputs x:{} ev:{} outputs y:{}  c:{}'.format(x, ev, f, c))
    return f, c, dict()

truemin = 0.39788735772973816