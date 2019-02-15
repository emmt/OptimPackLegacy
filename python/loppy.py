#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 10:30:00 2019

Use Cython optimpacklegacy to create a user friendly optimizer class

@author: rfetick
"""

import optimpacklegacy
import numpy as np

STATUS_OP_ERROR = -2
STATUS_NOT_FINISHED = -1
STATUS_WAIT = 0
STATUS_ITMAX = 1
STATUS_CONV = 2

TASK_START = 0
TASK_FG = 1
TASK_ISFREE = 2
TASK_NEWX = 3
TASK_CONV = 4
TASK_WARN = 5
TASK_ERROR = 6

class Optimizer(object):
    def __init__(self,func,grad,x,callback=None,blow=None,bup=None,verbose=False,
                 itmax=1000,fatol=1e-13,frtol=1e-11,**kwargs):
        self.func = func
        self.grad = grad
        self.kwargs = kwargs
        self.callback = callback
        
        self.x = np.copy(x) # Copy to avoid overwriting user input `x`
        self.n = len(x)
        self.m = 5
        
        size = optimpacklegacy.py_vmlmb_workspace_size(self.n,self.m)
        buffer = np.zeros(size,dtype=np.int8)
        self.ws = optimpacklegacy.py_vmlmb_monolithic_workspace_init(buffer,self.n,self.m)
        
        self.blow = blow
        self.bup = bup
        self.verbose = verbose
        self.itmax = int(itmax)
        
        self.fatol=np.float64(fatol)
        self.frtol=np.float64(frtol)
        
        optimpacklegacy.py_vmlmb_set_fatol(self.ws,fatol)
        optimpacklegacy.py_vmlmb_set_frtol(self.ws,frtol)
        
        self.last_f = None
        self.last_g = None
        self.task = 0
        self.iter = 0
        self.message = "(empty message)"
        self.status = STATUS_WAIT
        
        if (self.blow is not None) and (self.bup is not None):
            if np.any(self.blow>self.bup):
                raise ValueError("Some lower bounds are higher than upper bounds")

        # Eventually apply bounds on `x`
        self.applybounds()
		
        # Check types, last_f and last_g are populated for `active`
        self.checktypes()
        
        # Create `active` set, using last_g previously computed
        self.active = np.zeros(self.n,dtype=np.int32)
        self.activate()
        
        

    def __repr__(self):
        return "Optimpack optimizer [iter=%u] [task=%u] [status=%u]"%(self.iter,self.task,self.status)
    
    def activate(self):
        """Compute `active` set in case of bounds"""
        self.active.fill(1)
        if self.blow is not None:
            self.active *= (1-(self.x<=self.blow)*(self.last_g>=0)).astype(np.int32)
        if self.bup is not None:
            self.active *= (1-(self.x>=self.bup)*(self.last_g<=0)).astype(np.int32)
	
    def applybounds(self):
        """Apply bounds on `x` (projection step)"""
        if self.blow is not None:
            self.x = self.x*(self.x>self.blow) + self.blow*(self.x<=self.blow)
        if self.bup is not None:
            self.x = self.x*(self.x<self.bup) + self.bup*(self.x>=self.bup)

    def step(self):
        if self.task == TASK_START or self.task == TASK_FG:
            self.applybounds()
            self.last_f = self.func(self.x,**self.kwargs)
            self.last_g = self.grad(self.x,**self.kwargs)
        elif self.task == TASK_ISFREE:
            self.activate()
        elif self.task == TASK_NEWX:
            if self.callback is not None:
                self.callback(self.x)
        elif self.task == TASK_CONV:
            self.status = STATUS_CONV
        elif self.task == TASK_ERROR or self.task == TASK_WARN:
            self.status = STATUS_OP_ERROR
        
        self.task = optimpacklegacy.py_vmlmb_iterate(self.ws,self.x,self.last_f,self.last_g,self.active)
        self.iter += 1
        
    def run(self):
        self.status = STATUS_NOT_FINISHED
        while (self.iter<self.itmax) and (self.task not in [TASK_CONV,TASK_ERROR,TASK_WARN]):
            self.step()
            
        if self.iter>=self.itmax:
            self.status = STATUS_ITMAX
            
        return self.status


    def checktypes(self):
        """ Check if user inputs are okay"""
        try:
            f = self.func(self.x,**self.kwargs)
        except:
            print("Couldn't evaluate `func(x,**kwargs)`")
		
        try:
            g = self.grad(self.x,**self.kwargs)
        except:
            print("Couldn't evaluate `grad(x,**kwargs)`")
		
        self.last_f = f #save func
        self.last_g = g #save grad
		
        if not type(self.x) is np.ndarray:
            raise TypeError("`x` must be a numpy.ndarray")
		
        if not type(f) is np.ndarray:
            raise TypeError("`func` must return a numpy.ndarray of length=1")
		
        if not type(g) is np.ndarray:
            raise TypeError("`grad` must return a numpy.ndarray")
		
        if self.x.ndim != 1:
            raise TypeError("`x` must be an array of dimension ndim=1")
		
        if g.ndim != 1:
            raise TypeError("`grad` must return an array of dimension ndim=1")
		
        if len(f) != 1:
            raise TypeError("`func` must return a numpy.ndarray of length=1")
		
        if not self.x.dtype == np.float64:
            raise TypeError("`x` must be of dtype=numpy.float64")
		
        if not f.dtype == np.float64:
            raise TypeError("`func` must return an array of dtype=numpy.float64")
		
        if not g.dtype == np.float64:
            raise TypeError("`grad` must return an array of dtype=numpy.float64")
		
        if len(g) != len(self.x):
            raise ValueError("`grad` must return an array of same size as `x`")

