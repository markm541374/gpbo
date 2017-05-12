import numpy as np
from matplotlib import pyplot as plt
import gpbo

from gpbo.core.optutils import ballradsearch as brs

f,a = plt.subplots(1)



def pd(x):
    a.plot(x[0],x[1],'r.')
    theta = abs(np.arctan(x[0]/x[1]))
    if np.sqrt(x[0]**2+x[1]**2)>0.4*theta+0.025:
        return 0
    else:
        return 1

print brs(2,1.,pd,ndirs=200, lineSh=1e-6)
plt.savefig('/home/mark/Desktop/tmp.png')