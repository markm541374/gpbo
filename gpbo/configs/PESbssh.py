# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

#random search on a test function
import acquisitions
import reccomenders
import objectives
import optimize

import scipy as sp
import os

path = os.path.join(os.path.expanduser('~'),'Dropbox/workspace/GPshared/results/PESbssh')
aqfn,aqpara = acquisitions.PESbs
aqpara['lb']=[-1.,-1.]
aqpara['ub']=[1.,1.]
aqpara['ev']['s']=1e-12
cfn = objectives.cfaexp(1.,0.75)

aqpara['traincfn'] = True

stoppara= {'nmax':100}
stopfn = optimize.nstopfn

reccfn,reccpara = reccomenders.gpashin
reccpara['lb']=aqpara['lb']
reccpara['ub']=aqpara['ub']
reccpara['everyn']=1
reccpara['onlyafter']=aqpara['nrandinit']+1



ojfw,xmin,ymin = objectives.genbiasedmat52ojf(len(aqpara['lb']),aqpara['lb'],aqpara['ub'],0.5)
ojf = objectives.costfnwrap(ojfw,cfn)
if True:
    from matplotlib import pyplot as plt
    import time
    n = 100
    x_ = sp.linspace(-1,1,n)
    y_ = sp.linspace(-1,1,n)
    z_ = sp.empty([n,n])
    s_ = sp.empty([n,n])
    for i in xrange(n):
        for j in xrange(n):
            m_= ojfw(sp.array([y_[j],x_[i]]),**{'s':1e-9,'xa':0.,'d':[sp.NaN]})
            z_[i,j] = m_[0]

    fig, ax = plt.subplots( nrows=1, ncols=1 ,figsize=(10,10))
    CS = ax.contour(x_,y_,z_,30)
    ax.clabel(CS, inline=1, fontsize=8)
    
    ax.plot(xmin[0],xmin[1],'ro')
    fig.savefig(os.path.join(os.path.expanduser('~'),'Dropbox/debugoutput','truegeneratedobjective'+time.strftime('%d_%m_%y_%H:%M:%S')+'.png'))

ojfchar = {'dx':len(aqpara['lb']),'dev':len(aqpara['ev'])}