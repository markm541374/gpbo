import pickle
import scipy as sp
from scipy import linalg as spl
from scipy import stats as sps
from matplotlib import pyplot as plt
from matplotlib import patches
n=47


R,xmin = pickle.load(open("dbout/{}.p".format(n), "rb"))

f,ax = plt.subplots(2)


from sklearn import mixture
lowest_bic = sp.inf
bic = []
n_components_range = range(1, 7)
cv_types = ['spherical', 'tied', 'diag', 'full']
for cv_type in cv_types:
    for n_components in n_components_range:
        # Fit a Gaussian mixture with EM
        gmm = mixture.GaussianMixture(n_components=n_components,
                                          covariance_type=cv_type)
        gmm.fit(R)
        bic.append(gmm.bic(R))
        if bic[-1] < lowest_bic:
            lowest_bic = bic[-1]
            typ = cv_type
            best_gmm = gmm

clf = best_gmm
D_ = sp.inf
local = -1
for i,m in enumerate(clf.means_):
    dist=spl.norm(m-xmin)
    if dist<D_:
        local=i
        D_=dist

splot=ax[0]

Y_ = clf.predict(R)
for i, (mean, cov) in enumerate(zip(clf.means_, clf.covariances_)):
    color='blue' if i != local else 'red'
    if typ == 'spherical':
        cov = sp.diag([cov] * R.shape[1])
    elif typ == 'tied' or typ == 'diag':
        cov = sp.diag(cov)
    else:
        pass
    v, w = spl.eigh(cov)
    if not sp.any(Y_ == i):
        continue
    splot.scatter(R[Y_ == i, 0], R[Y_ == i, 1], .8, color=color)

    # Plot an ellipse to show the Gaussian component
    angle = sp.arctan2(w[0][1], w[0][0])
    angle = 180. * angle / sp.pi  # convert to degrees
    v = 2. * sp.sqrt(2.) * sp.sqrt(v)
    ell = patches.Ellipse(mean, v[0], v[1], 180. + angle, color=None, edgecolor=color,facecolor='none')
    splot.add_artist(ell)


#--------------------------------------------------
from sklearn import mixture
lowest_bic = sp.inf
bic = []
n_components_range = range(1, 7)
cv_types = ['spherical', 'tied', 'diag', 'full']
for cv_type in cv_types:
    for n_components in n_components_range:
        # Fit a Gaussian mixture with EM
        gmm = mixture.GaussianMixture(n_components=n_components,
                                          covariance_type=cv_type)
        gmm.fit(R)
        bic.append(gmm.bic(R))
        if bic[-1] < lowest_bic:
            lowest_bic = bic[-1]
            typ = cv_type
            best_gmm = gmm

clf = best_gmm
D_ = sp.inf
local = -1
print 'clf.means{}'.format(clf.means_)
for i,m in enumerate(clf.means_):
    dist=spl.norm(m-xmin)
    if dist<D_:
        local=i
        D_=dist
print 'local{}'.format(local)
locals = [local]
closestcov=clf.covariances_[local]
closestmean=clf.means_[local]
for i, (mean, cov) in enumerate(zip(clf.means_, clf.covariances_)):
    if typ == 'spherical':
        cov = sp.diag([cov] * R.shape[1])
    elif typ == 'tied' or typ == 'diag':
        cov = sp.diag(cov)
    else:
        pass
    print cov

    if (mean - closestmean).max() < sp.sqrt(cov.max()) and (mean - closestmean).max() < sp.sqrt(closestcov.max()):
        if closestcov.max() < 1.5 * cov.max() and closestcov.max() > 0.5 * cov.max():
            locals.append(i)
            print closestcov.max()
            print cov.max()
            print '-XXXX'
            print cov
            print closestcov
            print 'XXXX-'
print locals
splot=ax[1]

Y_ = clf.predict(R)
for i, (mean, cov) in enumerate(zip(clf.means_, clf.covariances_)):
    color='blue' if i not in locals else 'red'
    if typ == 'spherical':
        cov = sp.diag([cov] * R.shape[1])
    elif typ == 'tied' or typ == 'diag':
        cov = sp.diag(cov)
    else:
        pass
    v, w = spl.eigh(cov)
    if not sp.any(Y_ == i):
        continue
    splot.scatter(R[Y_ == i, 0], R[Y_ == i, 1], .8, color=color)

    # Plot an ellipse to show the Gaussian component
    angle = sp.arctan2(w[0][1], w[0][0])
    angle = 180. * angle / sp.pi  # convert to degrees
    v = 2. * sp.sqrt(2.) * sp.sqrt(v)
    ell = patches.Ellipse(mean, v[0], v[1], 180. + angle, color=None, edgecolor=color,facecolor='none')
    splot.add_artist(ell)




plt.show()

