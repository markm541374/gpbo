from __future__ import print_function
import scipy as sp



def f(x, **ev):
    #branin with a linear offset agains quadratic cost
    u = x[0] * 7.5 + 2.5
    v = x[1] * 7.5 + 2.5

    f = (-1.275 * (u / sp.pi) ** 2 + 5 * u / sp.pi + v - 6) ** 2 + (10. - 5. / (4 * sp.pi)) * sp.cos(u)+10.

    n = sp.random.normal() * sp.sqrt(ev['s'])
    if 'cheattrue' in ev.keys():
        if ev['cheattrue']:
            n = 0
    f+=n
    #let std =1e-2 > var = 1e-4 take 1+1 second. std=1e-4 will be 10000+1 ~3 hours
    c = (1e-4)/(ev['s']) + 1.

    print( 'f inputs x:{} ev:{} outputs y:{} (ytrue:{})  c:{}'.format(x, ev, f,f-n, c))
    return f, c, dict()

truemin = 0.39788735772973816
