import pickle
import scipy as sp
from scipy import linalg as spl
from scipy import stats as sps
from matplotlib import pyplot as plt
from matplotlib import patches
import gpbo
n=39


ER,M,V,Z_,Y_,R,Y,ymin,persist,Mf,Vf = pickle.load(open("dbout/{}.p".format(n), "rb"))
ns = Z_.shape[0]

P = sp.mean(Vf,axis=0)+sp.var(Mf,axis=0).reshape([1,ns])
V2 = Vf.max(axis=0)+sp.var(Mf,axis=0).reshape([1,ns])
M2 = Mf.min(axis=0)

ER2 = sp.empty([1,ns])
for i in xrange(ns):
    ER2[0,i] = gpbo.core.GPdc.EI(ymin,M2[i],V2[0,i])



count=[0,0,0]
zcount=[0,0,0]
for i in xrange(Z_.size):
    count[Z_[i]]+=1
    if ER[0,i]==0.:
        zcount[Z_[i]] += 1
        print '> >           {} | {} {} | {} {}'.format(ER2[0,i], Z_[i],Y[i, 1]-Y[i,0],(M[0,i]-ymin),(M[0,i]-ymin)/sp.sqrt(V[0,i]))
    elif Z_[i]==2:
        print '{} | {} {} | {} {}'.format(ER2[0,i], Z_[i],Y[i, 1]-Y[i,0],(M[0,i]-ymin),(M[0,i]-ymin)/sp.sqrt(V[0,i]))

print count
print zcount
print sp.sqrt(V[0,-1])

#print sp.hstack([Vf[:,:4].T,V[0,:4].reshape([4,1])])



