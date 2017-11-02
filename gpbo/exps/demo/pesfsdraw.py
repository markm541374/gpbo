import gpbo
import scipy as sp

D=4
s=1e-6
f,xmin,ymin = gpbo.core.objectives.genmat52ojf(D,[-1.]*D,[1.]*D,1.,[0.5]*D)
C=gpbo.core.config.pesfsdefault(f,D,80,s,'results','pesfs4d.csv')
out = gpbo.search(C)
print(out)