import scipy as sp
from gpbo.core import GPdc
from matplotlib import pyplot as plt
from gpbo.core.GPdc import kernel

X = sp.array([[-0.8],[-0.25],[0.25],[0.8]])
Y = sp.array([[0.3],[-0.2],[2.5],[5.]])
S = sp.array([[1e-3],[1e-3],[1e-3],[1e-3]])
D = [[sp.NaN],[sp.NaN],[sp.NaN],[sp.NaN]]
k= kernel(GPdc.SQUEXP, 1, sp.array([0.5, 0.2]))
g = GPdc.GPcore(X, Y, S, D,k)
ns=100
xaxis = sp.linspace(-1,1,ns)
m,v = g.infer_diag(xaxis,[[sp.NaN]]*ns)
s = sp.sqrt(v)
print (m-2.*s).flatten().shape
f,a = plt.subplots(2)
a[0].plot(X,Y,'r.')
a[0].fill_between(xaxis,(m-2.*s).flatten(),(m+2.*s).flatten(),facecolor='lightblue',edgecolor='lightblue')
a[0].plot(xaxis,m.flatten(),'b')

plt.show()
