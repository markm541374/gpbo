#!/usr/bin/env python2
#encoding: UTF-8

# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.


import optimize
import acquisitions
import reccomenders
import objectives

import scipy as sp
import os
import logging
logging.basicConfig(level=logging.DEBUG)


try:
    os.mkdir('test')
except OSError:
    pass
path = os.path.join(os.getcwd(),'test')

aqfn,aqpara = acquisitions.PESbs
#aqfn,aqpara = acquisitions.EIMAP
aqpara['lb']=[-1.,-1.]
aqpara['ub']=[1.,1.]
aqpara['ev']['s']=1e-12
cfn = objectives.cfaexp(1.,0.2)
ojf,xmin,ymin = objectives.genbiasedmat52ojf(len(aqpara['lb']),aqpara['lb'],aqpara['ub'],1.)
ojfn = objectives.costfnwrap(ojf,cfn)
ojfchar = {'dx':len(aqpara['lb']),'dev':len(aqpara['ev'])}
aqpara['cfn'] = cfn
aqpara['xau'] = 1.
aqpara['xal'] = 0.


#cfn = objectives.cf42
#ojfn,xmin,ymin = objectives.genmat52ojf(len(aqpara['lb']),aqpara['lb'],aqpara['ub'])
#ojfchar = {'dx':len(aqpara['lb']),'dev':len(aqpara['ev'])}



stoppara= {'nmax':20}
stopfn = optimize.nstopfn


reccfn,reccpara = reccomenders.gpasmap
reccpara['lb']=aqpara['lb']
reccpara['ub']=aqpara['ub']



O = optimize.optimizer(path,aqpara,aqfn,stoppara,stopfn,reccpara,reccfn,ojfn,ojfchar,checkrecc=True)

O.run()

