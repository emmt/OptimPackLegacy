#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 11:27:28 2019

@author: rfetick
"""

from loppy import Optimizer, decorate_fg
import numpy as np

#%% Define functions, start point and bounds
def banana(x,**kwargs):
    return 100*(x[1,...]-x[0,...]**2)**2 + (1.0-x[0,...])**2

def grad(x,**kwargs):
    u = x[1] - x[0]**2
    v = 1.0 - x[0]
    return np.array([-400.0*u*x[0]-2.0*v,200.0*u],dtype=np.float64)

x = np.zeros(2,dtype=np.float64)

bup = np.ones(2, dtype=np.float64)*.8
blow = np.zeros(2,dtype=np.float64)

#%% Run Optimizer

banana_and_grad = decorate_fg(banana,grad)

def callback(x,f,g):
    print("x = [%.3f,%.3f]    f = %.3f"%(x[0],x[1],f[0]))

lop = Optimizer(banana_and_grad,x,blow=blow,bup=bup,callback=callback)
lop.run()
