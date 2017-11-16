import pandas as pd
import numpy as np
import scipy as sp
from collections import defaultdict
def readset(files):
    D = defaultdict(lambda :[])
    for f in files:
        df = pd.read_csv(f)
        D['Ur'] = np.mean(np.ma.masked_invalid((1000-df['Um'].values-1)/df['Ut'].values))
        D['EIr'] = np.mean(np.ma.masked_invalid((1000-df['EIm'].values-1)/df['EIt'].values))
        D['CBr'] = np.mean(np.ma.masked_invalid((1000-df['CBm'].values-1)/df['CBt'].values))
        D['LRr'] = np.mean(np.ma.masked_invalid((1000-df['LRm'].values-1)/df['LRt'].values))
        D['LHr'] = np.mean(np.ma.masked_invalid((1000-df['LHm'].values-1)/df['LHt'].values))
        for k in df.keys():
            D[k].append(np.mean(np.ma.masked_invalid(df[k].values)))
    for k in D.keys():
        D[k] = np.mean(D[k])
    return D

def lines(D):
    kl = [D['Ukl'],D['EIkl'],D['CBkl'],D['LRkl'],D['LHkl']]
    ey = [D['Uy'],D['EIy'],D['CBy'],D['LRy'],D['LHy']]
    up = [(1000-f)/10. for f in D['Um'],D['EIm'],D['CBm'],D['LRm'],D['LHm']]
    ti = [D['Ut'],D['EIt'],D['CBt'],D['LRt'],D['LHt']]
    rp = [D['Ur'],D['EIr'],D['CBr'],D['LRr'],D['LHr']]
    return kl,ey,up,ti,rp

B2 = lines(readset(['results/b2support_{}.csv'.format(i) for i in range(5)]))
H3 = lines(readset(['results/h3support_{}.csv'.format(i) for i in range(5)]))
D4 = lines(readset(['results/d4support_{}.csv'.format(i) for i in range(5)]))
H6 = lines(readset(['results/h6support_{}.csv'.format(i) for i in range(5)]))

s = "branin 2d   &{:.3g} &{:.3g} &{:.3g} &{:.3g} &{:.3g} \\\\\n".format(*B2[0]) + \
    "hartmann 3d &{:.3g} &{:.3g} &{:.3g} &{:.3g} &{:.3g} \\\\\n".format(*H3[0]) + \
    "draw 4d     &{:.3g} &{:.3g} &{:.3g} &{:.3g} &{:.3g} \\\\\n".format(*D4[0]) + \
    "hartmann 6d &{:.3g} &{:.3g} &{:.3g} &{:.3g} &{:.3g} \\\\\n".format(*H6[0])
print(s+'\n')
s = "branin 2d   &{:.3g} &{:.3g} &{:.3g} &{:.3g} &{:.3g} \\\\\n".format(*B2[1]) + \
    "hartmann 3d &{:.3g} &{:.3g} &{:.3g} &{:.3g} &{:.3g} \\\\\n".format(*H3[1]) + \
    "draw 4d     &{:.3g} &{:.3g} &{:.3g} &{:.3g} &{:.3g} \\\\\n".format(*D4[1]) + \
    "hartmann 6d &{:.3g} &{:.3g} &{:.3g} &{:.3g} &{:.3g} \\\\\n".format(*H6[1])
print(s+'\n')
s = "branin 2d   &{:.1f} &{:.1f} &{:.1f} &{:.1f} &{:.1f} \\\\\n".format(*B2[2]) + \
    "hartmann 3d &{:.1f} &{:.1f} &{:.1f} &{:.1f} &{:.1f} \\\\\n".format(*H3[2]) + \
    "draw 4d     &{:.1f} &{:.1f} &{:.1f} &{:.1f} &{:.1f} \\\\\n".format(*D4[2]) + \
    "hartmann 6d &{:.1f} &{:.1f} &{:.1f} &{:.1f} &{:.1f} \\\\\n".format(*H6[2])
print(s+'\n')
s = "branin 2d   &{:.3g} &{:.3g} &{:.3g} &{:.3g} &{:.3g} \\\\\n".format(*B2[3]) + \
    "hartmann 3d &{:.3g} &{:.3g} &{:.3g} &{:.3g} &{:.3g} \\\\\n".format(*H3[3]) + \
    "draw 4d     &{:.3g} &{:.3g} &{:.3g} &{:.3g} &{:.3g} \\\\\n".format(*D4[3]) + \
    "hartmann 6d &{:.3g} &{:.3g} &{:.3g} &{:.3g} &{:.3g} \\\\\n".format(*H6[3])
print(s+'\n')
s = "branin 2d   &{:.3g} &{:.3g} &{:.3g} &{:.3g} &{:.3g} \\\\\n".format(*B2[4]) + \
    "hartmann 3d &{:.3g} &{:.3g} &{:.3g} &{:.3g} &{:.3g} \\\\\n".format(*H3[4]) + \
    "draw 4d     &{:.3g} &{:.3g} &{:.3g} &{:.3g} &{:.3g} \\\\\n".format(*D4[4]) + \
    "hartmann 6d &{:.3g} &{:.3g} &{:.3g} &{:.3g} &{:.3g} \\\\\n".format(*H6[4])
print(s+'\n')
