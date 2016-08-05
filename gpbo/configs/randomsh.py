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

path = os.path.join(os.path.expanduser('~'),'Dropbox/workspace/GPshared/results/randomsh')
aqfn,aqpara = acquisitions.random
aqpara['lb']=[-1.,-1.]
aqpara['ub']=[1.,1.]

stoppara= {'nmax':50}
stopfn = optimize.nstopfn

reccfn,reccpara = reccomenders.gpmap
reccpara['lb']=aqpara['lb']
reccpara['ub']=aqpara['ub']

cfn = objectives.cf42

ojf = objectives.trivialojf
ojfchar = {'dx':2,'dev':len(aqpara['ev'])}