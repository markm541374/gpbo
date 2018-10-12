import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
import copy
def kernel(x0,y0,x,y,L):
	return np.exp(-0.5*((x-x0)**2 + (y-y0)**2)/L**2)
def baseplot():
	f,a = plt.subplots(1,figsize=[5,5])
	n=160
	X = np.linspace(0,1,n)
	M = np.zeros([n,n])
	for i in range(n):
		for j in range(n):
			W = np.array([[0.27,0.26,0.25,0.24,0.23,0.22,0.21,0.28,0.35],
                                      [0.2,0.3,0.4,0.5,0.6,0.67,0.75,0.7,0.7],
                                      [1.,1.,1.,1.25,1.,1.,1.,1.,1.]])
			for p in range(W.shape[1]):
				M[i,j] += W[2,p]*kernel(W[1,p],W[0,p],X[i],0.5*X[j]+0.1,0.055)
	a.contour(X,X,0.0001*np.log(M)+M)
	a.set_xticklabels([])
	a.set_yticklabels([])
	return f,a

def grid(f,a):
	N = 8
	for i in range(N):
		for j in range(N):
			a.plot((i+0.5)/float(N),(j+0.5)/float(N),'ro')
	f.suptitle("Regular Grid",y=.995)
	return f,a

def slice(f,a):
	x0,y0 = 0.22,0.59
	a.plot(x0,y0,'ro')
	xoff = np.array([0.04,0.08,0.16,0.32,0.24,0.20,-0.04,-0.08,-0.06])
	a.plot(x0+xoff,y0*np.ones_like(xoff),'go')
	x1,y1=0.35,y0
	yoff = np.array([0.,-0.04,-0.08,-0.16,-0.32,-0.48,-0.4,0.04,0.08,0.16,0.32,0.24,0.20])
	a.plot(x1*np.ones_like(yoff),y0+yoff,'go')
	x2,y2 = x1,y1-0.2
	a.plot(x2,y2,'ro')
	a.arrow(x0,y0,x2-x0,y2-y0,width=0.005,head_width=0.03,length_includes_head=True,color='k')
	f.suptitle("Slice Sampling",y=0.995)
	return f,a
def square(x0,x1,y0,y1,a):
	t=0.005
	c='k'
	a.plot([x0,x1,x1,x0,x0],[y0,y0,y1,y1,y0],color=c,lw=t)
	xm,ym = 0.5*(x1+x0),0.5*(y0+y1)
	a.plot([x0,x0,x0,xm,xm,xm,x1,x1,x1],[y0,ym,y1,y0,ym,y1,y0,ym,y1],'ro')
	return

def quad(f,a):
	square(0,1,0,0.5,a)
	square(0,1,0.5,1,a)
	
	square(0,0.5,0,0.5,a)
	square(0.5,1,0,0.5,a)

	square(0.5,1.,0.5,1,a)
	square(0,0.5,0.5,1,a)

	square(0,0.5,0.75,1,a)
	square(0,0.5,0.5,0.75,a)

	square(0.25,0.5,0.5,0.75,a)
	square(0.,0.25,0.5,0.75,a)

	square(0.25,0.5,0.5,0.625,a)
	square(0.25,0.5,0.625,0.75,a)

	square(0.25,0.375,0.625,0.75,a)
	square(0.375,0.5,0.625,0.75,a)

	square(0.3125,0.375,0.625,0.75,a)
	square(0.25,0.3125,0.625,0.75,a)
	f.suptitle("Adaptive Quadrature",y=0.995)
	return f,a
if __name__=="__main__":
	#f,a = baseplot()
	#f.savefig("base.pdf")
	
	f1,a1 = grid(*baseplot())
	plt.tight_layout()
	f1.savefig("figs/grid.pdf")
	f2,a2 = slice(*baseplot())
	plt.tight_layout()
	f2.savefig("figs/slice.pdf")
	f3,a3 = quad(*baseplot())
	plt.tight_layout(rect=[0,0,1,1])
	f3.savefig("figs/quad.pdf")
	
	
