#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 11:27:28 2019

@author: rfetick
"""

from loppy import Optimizer
import numpy as np

#%% Define functions, start point and bounds

def banana_array(x,**kwargs):
    return 100*(x[1,...]-x[0,...]**2)**2 + (1.0-x[0,...])**2

def banana(x,**kwargs):
    a = np.zeros(1,dtype=np.float64)
    a[0] = banana_array(x,**kwargs)
    return a

def banana_grad(x,**kwargs):
    u = x[1] - x[0]**2
    v = 1.0 - x[0]
    return np.array([-400.0*u*x[0]-2.0*v,200.0*u],dtype=np.float64)

x = np.zeros(2,dtype=np.float64)

bup = np.ones(2, dtype=np.float64)*.8
blow = np.zeros(2,dtype=np.float64)

#%% Run Optimizer

def callback(x):
    print("x = [%.3f,%.3f]"%(x[0],x[1]))

lop = Optimizer(banana,banana_grad,x,blow=blow,bup=bup,callback=callback)
lop.run()
