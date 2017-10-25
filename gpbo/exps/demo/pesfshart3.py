import gpbo
import scipy as sp


#gpbo.core.debugoutput['datavis']=True


D=3
n=50
s=1e-6
f = gpbo.core.objectives.shifthart3

C=gpbo.core.config.pesfsdefault(f,D,150,s,'results','pesfs.csv')
#C.aqfn = gpbo.core.acquisitions.vmaxaq
out = gpbo.search(C)
print(out)