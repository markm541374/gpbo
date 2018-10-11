import scipy as sp
from scipy import stats
import numpy as np
from matplotlib import pyplot as plt
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
import tqdm

n=10000
m=20
idx = np.linspace(5,n,3,dtype=np.int)
x = np.linspace(-4,4,201)
#X0 = np.random.normal(loc=0.,scale=1.,size=n)
J0 = np.empty(m)
for i in tqdm.tqdm(range(m)):
    X0 = np.random.normal(loc=0.,scale=1.,size=n)
    J0[i]=sp.stats.jarque_bera(X0)[0]
#JS0 = np.empty([idx.size,2])
#for k in tqdm.tqdm(range(idx.size)):
#    JS0[k,:] = sp.stats.jarque_bera(X0[:idx[k]])
p0 = sp.stats.norm.pdf(x,loc=0.,scale=1.)

J1 = np.empty(m)
for i in tqdm.tqdm(range(m)):
    S1 = np.vstack([np.random.normal(loc=0.2,scale=1.1,size=n),np.random.normal(loc=-0.2,scale=0.9,size=n)])
    I = sp.random.randint(2,size=n)
    X1 = np.array([S1[I[k],k] for k in range(S1.shape[1])])
    J1[i]=sp.stats.jarque_bera(X1)[0]
#S1 = np.vstack([np.random.normal(loc=0.1,scale=1.1,size=n),np.random.normal(loc=-0.1,scale=0.9,size=n)])
#X1 = S1[sp.random.randint(2,size=n),:]
#JS1 = np.empty([idx.size,2])
#for k in tqdm.tqdm(range(idx.size)):
#    JS1[k,:] = sp.stats.jarque_bera(X1[:idx[k]])
#J1 = sp.stats.jarque_bera(X1)
p1 = 0.5*sp.stats.norm.pdf(x,loc=0.2,scale=1.1)+0.5*sp.stats.norm.pdf(x,loc=-0.2,scale=0.9)

print(J0, np.mean(J0))
print(J1,np.mean(J1))

#S2 = np.vstack([np.random.normal(loc=0.9,scale=2.1,size=n),np.random.normal(loc=-0.1,scale=0.9,size=n)])
#X2 = S2[sp.random.randint(2,size=n),:]
#JS2 = np.empty([idx.size,2])
#for k in tqdm.tqdm(range(idx.size)):
#    JS2[k,:] = sp.stats.jarque_bera(X2[:idx[k]])
#J2 = sp.stats.jarque_bera(X2)

#p2 = 0.5*sp.stats.norm.pdf(x,loc=0.9,scale=2.1)+0.5*sp.stats.norm.pdf(x,loc=-0.1,scale=0.9)
"""
f,a = plt.subplots(1)
at = a.twinx()
a.plot(idx, JS0[:,0],'b')
at.plot(idx, JS0[:,1],'b.')
a.plot(idx, JS1[:,0],'r')
at.plot(idx, JS1[:,1],'r.')
a.plot(idx, JS2[:,0],'g')
at.plot(idx, JS2[:,1],'g.')
a.set_yscale('log')
"""

f2,a2 = plt.subplots()
a2.plot(x,p0,colors[1],label="Unit Normal")
a2.plot(x,p1,colors[2],label="Combination of Normal")
a2.set_xlabel("$x$")
a2.set_ylabel("$p(x)$")
a2.legend(loc=2)
a2.set_title('Distribution Comparison')
f2.savefig('figs/refdist.pdf')

plt.show()
