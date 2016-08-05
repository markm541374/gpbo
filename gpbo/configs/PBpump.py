# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

#random search on a test function
import gpbo.core.acquisitions as acquisitions
import gpbo.core.reccomenders as reccomenders
import gpbo.core.objectives as objectives
import gpbo.core.optimize as optimize
import scipy as sp
import os

path = os.path.join(os.path.expanduser('~'),'Dropbox/workspace/GPshared/results/PBpump')
aqfn,aqpara = acquisitions.PESbs
aqpara['lb']=[-1.,-1.,-1.]
aqpara['ub']=[1.,1.,1.]
aqpara['ev']['s']=1e-6
aqpara['traincfn'] = True
aqpara['mprior']=sp.array([1.,0.,0.,0.,-2.])
aqpara['sprior']=sp.array([1.,1.,1.,1.,1.])

stoppara= {'nmax':100}
stopfn = optimize.nstopfn

reccfn,reccpara = reccomenders.gpasmap
reccpara['lb']=aqpara['lb']
reccpara['ub']=aqpara['ub']
reccpara['everyn']=1
reccpara['onlyafter']=aqpara['nrandinit']
reccpara['mprior'] = aqpara['mprior']
reccpara['sprior'] = aqpara['sprior']
import sys
sys.path.append(os.path.join(os.path.expanduser('~'),'Dropbox/workspace/control'))
import opt
ojf = opt.BOSojf
ojfchar = {'dx':len(aqpara['lb']),'dev':len(aqpara['ev'])}
