import gpbo
import scipy as sp

D=3
n=50
s=1e-6
f = gpbo.core.objectives.shifthart3

C=gpbo.core.config.pesfsdefault(f,D,120,s,'results','pesfsh.csv')
out = gpbo.search(C)
print(out)