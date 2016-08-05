# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.
def run():
    import gpbo.core.optimize
    import gpbo.core.acquisitions
    import gpbo.core.reccomenders
    import gpbo.core.objectives

    import scipy as sp
    import os
    import sys

    import logging
    logging.basicConfig(level=logging.DEBUG)

    sys.path.append('configs')
    #import gpbo.configs.randomsh as optconfig
    import gpbo.configs.gridsh as optconfig
    #import gpbo.configs.EIMLsh as optconfig
    #import gpbo.configs.PESbssh as optconfig
    #import gpbo.configs.EIpumppid as optconfig
    #import gpbo.configs.PESpump as optconfig
    #import gpbo.configs.PBpump as optconfig

    O = gpbo.core.optimize.optimizer(optconfig.path, optconfig.aqpara, optconfig.aqfn, optconfig.stoppara, optconfig.stopfn, optconfig.reccpara, optconfig.reccfn, optconfig.ojf, optconfig.ojfchar, checkrecc=True)

    O.run()
    return

if __name__=="__main__":
    run()