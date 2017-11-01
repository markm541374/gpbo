import gpbo
import scipy as sp
from matplotlib import pyplot as plt

#gpbo.core.debugoutput['datavis']=True


D=4
s=1e-6
f,xmin,ymin = gpbo.core.objectives.genmat52ojf(D,[-1.]*D,[1.]*D,1.,[0.5]*D)
"""
fig,a = plt.subplots()
l = 200
x_ = sp.linspace(-1,1,l)
y_ = sp.linspace(-1,1,l)
z_ = sp.empty([l,l])
s_ = sp.empty([l,l])
for i in range(l):
    for j in range(l):
        m_,_,__ = f([y_[j],x_[i],0.],**{'d':[sp.NaN],'s':1e-9})
        z_[i,j] = m_
CS = a.contour(x_,y_,z_,20)
a.clabel(CS, inline=1, fontsize=10)
fig.savefig('dbout/ojf.svg')
"""
C=gpbo.core.config.pesfsdefault(f,D,80,s,'results','pesfs4d.csv')
out = gpbo.search(C)
print(out)