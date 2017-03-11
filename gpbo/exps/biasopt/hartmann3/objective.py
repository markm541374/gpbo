from __future__ import print_function
import scipy as sp
import numpy as np


def f(x, **ev):
    #hartmann4 with a linear offset agains quadratic cost
    sn = ev['xa']
    shift = [0.15,-0.2,-0.3]
    z = [0.5*i+0.5+sn*shift[j] for j,i in enumerate(x)]

    al = np.array([1.,1.2,3.,3.2]).T
    A = np.array([[3.,   10., 30.],
                  [0.1,  10., 35.],
                  [3.,   10., 30.],
                  [0.1,  10., 35.]])
    P = 0.0001*np.array([[3689, 1170, 2673],
                         [4699, 4387, 7470],
                         [1091, 8732, 5547],
                         [381, 5743, 8828]])

    outer = 0
    for ii in range(4):
        inner = 0
        for jj in range(3):
            xj = z[jj]
            Aij = A[ii, jj]
            Pij = P[ii, jj]
            inner += Aij * (xj - Pij)**2



        new = al[ii] * np.exp(-inner)
        outer = outer + new
    f = -outer + sn*0.5
    c = 120. + 1680. * (1. - sn) ** 2

    print( 'f inputs x:{} ev:{} outputs y:{}  c:{}'.format(x, ev, f, c))
    return f, c, dict()

truemin = -3.86278