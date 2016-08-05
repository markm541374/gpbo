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

path = os.path.join(os.path.expanduser('~'),'Dropbox/workspace/GPshared/results/EIpump')
aqfn,aqpara = acquisitions.EIMAP
aqpara['lb']=[-2.,-2.,-2.]
aqpara['ub']=[1.,1.,1.]

stoppara= {'nmax':100}
stopfn = optimize.nstopfn

reccfn,reccpara = reccomenders.gpmap
reccpara['lb']=aqpara['lb']
reccpara['ub']=aqpara['ub']
reccpara['everyn']=4
reccpara['onlyafter']=aqpara['nrandinit']

cfn = objectives.cf42

import sys
sys.path.append(os.path.join(os.path.expanduser('~'),'Dropbox/workspace/control'))
import opt

ojf = opt.BOojf
ojfchar = {'dx':len(aqpara['lb']),'dev':len(aqpara['ev'])}