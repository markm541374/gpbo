import gpbo
import scipy as sp

D=6
n=50
s=1e-6


f = gpbo.core.objectives.shifthart6

C=gpbo.core.config.pesfsdefault(f,D,120,s,'results','pesfsh6.csv')
#C.aqfn = gpbo.core.acquisitions.vmaxaq
out = gpbo.search(C)
print(out)