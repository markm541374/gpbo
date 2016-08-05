#!/usr/bin/env python2
#encoding: UTF-8

# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.
import sys
import os
sys.path.append("./GPpython/src")
sys.path.append("./GPpython/src/tests")

if __name__ == "__main__":
    #GP demo
    import testdraw
    import testGPdc
    import testdrawhyp
    
    #compare PES with VF (red)  to PES (dots )and EI (x) at fixed steps on branin with an impoesd cost/noise relation
    
    os.chdir('GPpython/src/')
    #I know this is bad but it's easier than fixing the relative paths
    os.system('python  expbraninnoise.py')
    
    #compare PES with EV (blue) to PES and various fixed steps on synthetic functions. without the cache of results this will take a very long time. The smaller steps converge to a biased results as they are not in the s=0 plane
    os.system('python  expIPS.py')