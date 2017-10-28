import gpbo
import scipy as sp


#gpbo.core.debugoutput['datavis']=True


D=2
s=1e-6
f = gpbo.core.objectives.shiftbraninojf

C=gpbo.core.config.pesfsdefault(f,D,150,s,'results','pesfsb.csv')
#C.aqfn = gpbo.core.acquisitions.vmaxaq
out = gpbo.search(C)
print(out)