# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

import optimize
import acquisitions
import reccomenders
import objectives

import scipy as sp
import os
import sys

import logging
logging.basicConfig(level=logging.DEBUG)

sys.path.append('configs')
#import randomsh as optconfig
#import gridsh as optconfig
#import EIMLsh as optconfig
import PESbssh as optconfig
#import EIpumppid as optconfig
#import PESpump as optconfig
#import PBpump as optconfig

O = optimize.optimizer(optconfig.path,optconfig.aqpara,optconfig.aqfn,optconfig.stoppara,optconfig.stopfn,optconfig.reccpara,optconfig.reccfn,optconfig.ojf,optconfig.ojfchar,checkrecc=True)

O.run()
