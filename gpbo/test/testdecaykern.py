import scipy as sp
from gpbo.core import GPdc
from matplotlib import pyplot as plt
from gpbo.core.GPdc import kernel

xs=sp.linspace(0.,35,10)
X = sp.array([xs]).T
Y = sp.array([map(lambda x:sp.exp(-0.5*x)-1.+0.01*sp.random.normal(),xs)]).T
S = sp.array([[1e-2]*X.shape[0]])
D = [[sp.NaN]]*X.shape[0]
ns=100
xaxis = sp.linspace(0.,60.,ns)
f,a = plt.subplots(2)

#MAPHYP = GPdc.searchMAPhyp(X, Y, S, D, sp.array([-1., -1.]), sp.array([3., 3.]), GPdc.DEC1)
#k = kernel(GPdc.DEC1,1,MAPHYP)
#print MAPHYP
k = kernel(GPdc.CPDEC1,1,sp.array([100.,200.,1.]))
g = GPdc.GPcore(X, Y, S, D,k)
m,v = g.infer_diag(xaxis,[[sp.NaN]]*ns)
s = sp.sqrt(v)
a[0].plot(X,Y,'r.')
a[0].fill_between(xaxis,(m-2.*s).flatten(),(m+2.*s).flatten(),facecolor='lightblue',edgecolor='lightblue')
a[0].plot(xaxis,m.flatten(),'b')
#for i in xrange(nd):
#    a.plot(xaxis,R[i,:],'g')


#
plt.show()
