from __future__ import print_function
import scipy as sp
import numpy as np


def f(x, **ev):
    #hartmann4 with a linear offset agains quadratic cost
    sn = ev['xa']
    shift = [0.1,-0.15,0.2,0.05,-0.25,0.1]
    z = [0.5*i+0.5+sn*shift[j] for j,i in enumerate(x)]

    al = np.array([1.,1.2,3.,3.2]).T
    A = np.array([[10., 3., 17., 3.5, 1.7, 8.],
                  [0.05, 10., 17., 0.1, 8, 14.],
                  [3., 3.5, 1.7, 10., 17., 8.],
                  [17., 8., 0.05, 10., 0.1, 14.]])
    P = 0.0001 * np.array([[1312., 1696., 5569., 124., 8283., 5886.],
                           [2329., 4135., 8307., 3736., 1004., 9991.],
                           [2348., 1451., 3522., 2883., 3047., 6650.],
                           [4047., 8828., 8732., 5743., 1091., 381.]])

    outer = 0
    for ii in range(4):
        inner = 0
        for jj in range(6):
            xj = z[jj]
            Aij = A[ii, jj]
            Pij = P[ii, jj]
            inner += Aij * (xj - Pij)**2



        new = al[ii] * np.exp(-inner)
        outer = outer + new
    f = -outer + sn*0.5
    c = 30. + 570. * (1. - sn) ** 2

    print( 'f inputs x:{} ev:{} outputs y:{}  c:{}'.format(x, ev, f, c))
    return f, c, dict()

truemin = -3.32237
xtruemin = [0.20169,0.150011,0.476874,0.275332,0.311652,0.6573]

if __name__=="__main__":
    print(f([-0.598079561, -0.6995884774, -0.0438957476, -0.4499314129, -0.3758573388, 0.3127572016],**{'xa':0.}))

