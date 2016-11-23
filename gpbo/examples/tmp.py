import pickle
import scipy as sp
from scipy import linalg as spl
from scipy import stats as sps
from matplotlib import pyplot as plt
from matplotlib import patches
n=90


W,R,Y,ER,Y_,localset,CommitRegret,IRegret,ORegret,ERegret,n,persist,means, covs,M,V,xmin,ymin,W_,A = pickle.load(open("dbout/{}.p".format(n), "rb"))
fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(27, 40))

ax[0,1].plot(W[:, 0], W[:, 1], 'bx')
ax[0,1].plot(R[:, 0], R[:, 1], 'r.')

fullbound=ER.sum()
innerbound=0.
outerbound=0.
for i in xrange(ER.shape[1]):
    if i in A:
        j=A.index(i)
        if Y_[j] in localset:
            innerbound+=ER[0,i]
            print 'inAinY {}'.format(ER[0,i])
        else:
            outerbound+=ER[0,i]
            print 'inA outY ER {} DR {}'.format(ER[0,i],Y[j, 1] - Y[j, 0])
    else:
        outerbound+=ER[0,i]
        print 'outA {}'.format(ER[0, i])
        if ER[0,i]>0.:
            ax[1,0].plot(W[i,0],W[i,1],'.')

msg = 'n            : {}\n' \
      'localset     : {}\n\n' \
      'Com Regret   : {}\n' \
      'Com Regret B : {}\n\n' \
      'Loc Regret   : {}\n' \
      'Loc Regret B : {}\n\n' \
      'All Regret   : {}\n' \
      'All Regret B : {}\n' \
      ''.format(n,localset,CommitRegret,outerbound,IRegret,innerbound,ERegret,fullbound)
l_x, u_x = ax[1, 1].get_xlim()
l_y, u_y = ax[1, 1].get_ylim()
ax[1,1].text(l_x + 0.02 * (u_x - l_x), l_y + 0.02 * (u_y - l_y), '\n----\n'+msg,fontdict={'size':12})

for i, (mean, cov) in enumerate(zip(means, covs)):
    color = 'blue' if i not in localset else 'red'

    v, w = spl.eigh(cov)
    # Plot an ellipse to show the Gaussian component
    angle = sp.arctan2(w[0][1], w[0][0])
    angle = 180. * angle / sp.pi  # convert to degrees
    v = 2. * sp.sqrt(2.) * sp.sqrt(v)
    ell = patches.Ellipse(mean, v[0], v[1], 180. + angle, color=None, edgecolor=color, facecolor='none')
    ax[0,2].add_artist(ell)
    ax[0,2].text(mean[0],mean[1],str(i),fontdict={'color':color})
ax[0,2].axis([-1,1,-1,1])


plt.show()