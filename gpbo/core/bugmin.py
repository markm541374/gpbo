#!/usr/bin/env python2
#encoding: UTF-8

# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

import OPTutils
import scipy as sp
from matplotlib import pyplot as plt
import GPdc

import eprop
import pickle
[m0,V0,Y,Z,F,z] = pickle.load(open('epkill.p','rb'))
print m0
print V0
print Y
print Z
print F
print z

import eprop
eprop.expectation_prop(m0,V0,Y,Z,F,z)