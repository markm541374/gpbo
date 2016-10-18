
import scipy as sp
import gpbo

def testsimple(ns,no):
    import fab
    xl = [-2.,-2.]
    xu = [2.,2.]
    sworst=1.
    sbest=0.
    def f(x, s):
        y = -sp.cos(x[0]*0.5) -sp.cos(0.5*x[1])+0.1*s**2+2.
        c = 1-0.5*s
        return y, c
    for i in xrange(no):
        fab.runfab(xl,xu,sworst,sbest,f,ns,2,ofpath='results/fabout{}.csv'.format(i))
    return


run=False
plot=True

nopts=10


if run:
    testsimple(42,nopts)
if plot:
    from matplotlib import pyplot as plt

    d = [gpbo.optimize.readoptdata('results/fabout{}.csv'.format(k)) for k in xrange(nopts)]

    f,a=plt.subplots(2)
    for i in xrange(nopts):
        a[0].plot(d[i]['cacc'],d[i]['trueyatxrecc'],'b')
        a[1].plot(d[i]['cacc'],d[i]['c'])
    a[0].set_yscale('log')
    f.savefig('plots/out0.pdf')
    f,a=plt.subplots(1)
    xaxis = sp.linspace(0, 100, 100)

    low0, med0, upp0 = gpbo.core.ESutils.quartsirregular([d[k]['cacc'] for k in xrange(nopts)],
                                                         [d[k]['trueyatxrecc'] for k in xrange(nopts)], xaxis)

    a.fill_between(xaxis, low0, upp0, facecolor='lightblue', edgecolor='lightblue', alpha=0.5)
    a.plot(xaxis, med0, 'b')
    a.set_yscale('log')
    f.savefig('plots/out1.pdf')
    plt.show()
