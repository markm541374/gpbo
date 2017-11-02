import gpbo
import scipy as sp


D=2
s=1e-6
f = gpbo.core.objectives.shiftbraninojf

C=gpbo.core.config.pesfsdefault(f,D,150,s,'results','pesfsb.csv')
out = gpbo.search(C)
print(out)