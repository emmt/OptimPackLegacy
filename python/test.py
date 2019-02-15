#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 08 17:18:00 2019

Test optimpacklegacy in Python

@author: rfetick
"""

import optimpacklegacy
import numpy as np

## DEFINE FUNCTION TO MINIMZIE ##

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
isfree = np.zeros(2,dtype=np.int32)

bup = np.ones(2, dtype=np.float64)*.8
blow = np.zeros(2,dtype=np.float64)

## INIT OPTIMPACKLEGACY ##

n = len(x)
m = 5

size = optimpacklegacy.py_vmlmb_workspace_size(n,m)
buffer = np.zeros(size,dtype=np.int8)
ws = optimpacklegacy.py_vmlmb_monolithic_workspace_init(buffer,n,m)

## LOOP ##

task = 1
nb_eval = 0

for i in range(1000):
    if task == 1:
        # apply bounds on x
        x = x*(x>blow)+blow*(x<=blow)
        x = x*(x<bup)+bup*(x>=bup)
        
        f=banana(x)
        g=banana_grad(x)
        nb_eval += 1
    if task == 2:
        # compute free variables
        # (X[i] > XMIN[i] || G[i] < 0) && (X[i] < XMAX[i] || G[i] > 0)
        isfree = (((x>blow) + (g<0))*((x<bup)+(g>0))).astype(np.int32)
        
    if task == 3:
        print("x = [%.2f,%.2f]"%tuple(x))
    if task == 4:
        print("Algo converged in %u iter and %u evaluations"%(i,nb_eval))
        break
    if task>4:
        raise ValueError("Some error occured")
    task = optimpacklegacy.py_vmlmb_iterate(ws,x,f,g,isfree)
    


