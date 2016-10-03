
import gpbo.core.objectives as objectives
import gpbo.core.optimize as optimize

import scipy as sp
import os
import time

D=2
lb = [-1., -1.]
ub = [1., 1.]

cfn = objectives.cfaexp(1., 3.)

stoppara = {'nmax': 100}
stopfn = optimize.nstopfn


ojfw, xmin, ymin = objectives.gendecayingpositiveojf(D, lb, ub)
ojf = objectives.costfnwrap(ojfw, cfn)


#show the generated objective
if True:
    from matplotlib import pyplot as plt


    n = 100
    x_ = sp.linspace(-1, 1, n)
    y_ = sp.linspace(-1, 1, n)
    z_ = sp.empty([n, n])
    s_ = sp.empty([n, n])
    for i in xrange(n):
        for j in xrange(n):
            m_ = ojf(sp.array([y_[j], x_[i]]), **{'s': 0, 'xa': 0., 'd': [sp.NaN]})
            z_[i, j] = m_[0]

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))
    CS = ax.contour(x_, y_, z_, 30)
    ax.clabel(CS, inline=1, fontsize=8)

    ax.plot(xmin[0], xmin[1], 'ro')

    n = 100
    x_ = sp.linspace(-1, 1, n)
    y_ = sp.linspace(-1, 1, n)
    z_ = sp.empty([n, n])
    s_ = sp.empty([n, n])
    for i in xrange(n):
        for j in xrange(n):
            m_ = ojf(sp.array([y_[j], x_[i]]), **{'s': 0, 'xa': 1., 'd': [sp.NaN]})
            z_[i, j] = m_[0]

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))
    CS = ax.contour(x_, y_, z_, 30)
    ax.clabel(CS, inline=1, fontsize=8)

    n = 100
    x_ = sp.linspace(0, 1, n)
    y_ = sp.linspace(-1, 1, n)
    z_ = sp.empty([n, n])
    s_ = sp.empty([n, n])
    for i in xrange(n):
        for j in xrange(n):
            m_ = ojf(sp.array([0,y_[j]]), **{'s': 0, 'xa': x_[i], 'd': [sp.NaN]})
            z_[j, i] = m_[0]

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))
    CS = ax.contour(x_, y_, z_, 30)
    ax.clabel(CS, inline=1, fontsize=8)


    plt.show()
    #fig.savefig(os.path.join('.', 'truegeneratedobjective' + time.strftime('%d_%m_%y_%H:%M:%S') + '.png'))
    del (fig)
