import fab
import scipy as sp

def testsimple():
    xl = [-2.,-2.]
    xu = [2.,2.]
    sworst=1.
    sbest=0.
    def f(x, s):
        y = -sp.cos(x[0]*0.5) -sp.cos(0.5*x[1])+0.1*s**2+2.
        c = 1-0.5*s
        return y, c

    fab.runfab(xl,xu,sworst,sbest,f,44,2,ofpath='fabout.csv')
    return

if __name__=="__main__":
    testsimple()