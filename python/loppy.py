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


def decorate_array(func):
    """
    Decorator to use if the function `func` returns a scalar, and makes it return an array of size=1
    """
    def decorated(*args,**kwargs):
        out = np.zeros(1,dtype=np.float64)
        out[0] = func(*args,**kwargs)
        return out
    return decorated


class Optimizer(object):
    """
    Solves the following minimisation problem, with eventual bounds:
    	xs = argmin{f(x,**kwargs)} , b_low < x <b_up
    
    Notes
    -----
    The minimizer uses VMLM-B method (Eric Thiébaut - CRAL, Lyon, France) and requires computation of the gradient by the user. It is designed to solve problems with high number of unknowns.
    
    The core of the algorithm is based on the OptimPackLegacy library written in C, under GPL v2 license. Adaptation to Python was made using Cython
    
    Example
    -------
    >>> from loppy import Optimizer
    >>> func = ... # define your function
    >>> grad = ... # define its gradient
    >>> x0 = ...   # give an initial guess
    >>> opti = Optimizer(func,grad,x0,blow=...,bup=...)
    >>> opti.run()
    >>> print(opti.x)
    
    See also
    --------
    Optimizer.__init__.__doc__
    
    """
    def __init__(self,func,grad,x,callback=None,blow=None,bup=None,
                 verbose=False,itmax=1000,fatol=None,frtol=None,fmin=None,
                 delta=None,epsilon=None,**kwargs):
        """
    Optimizer initialisation method
    	
    Parameters
    ----------
    func : callable
        Function to minimize
    grad : callable
        Gradient of the function to minimize
    x : numpy.ndarray(dtype=float64)
        Initial guess on parameters
       	
    Keywords
    --------
    callback : callable
        Function to call after a new `x` is found
        (e.g. for visualization)
    blow : scalar or numpy.ndarray(dtype=numpy.float64)
        Lower bounds on 'x' (default: no bounds)
    bup : scalar or numpy.ndarray(dtype=numpy.float64)
        Upper bounds on 'x' (default: no bounds)
    verbose : bool
        Display inline information about minimization (default=False)
    itmax : int
        Number maximum of OptimPack iterations (default=1000)
    fmin : float
        Minimal expected value of `func`
    fatol : float
        Absolute tolerance on f such as
        |f(x_{k})-f(x_{k+1})| < fatol
    frtol : float
        Relative tolerance on f such as
        2*|f(x_{k})-f(x_{k+1})|/(f(x_{k})+f(x_{k+1})|) < frtol
    delta : float
    epsilon : float
    **kwargs : set of keywords
        Eventual keywords to be given to your functions `func` and `grad`
        """
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
        
        if fmin is not None:
            optimpacklegacy.py_vmlmb_set_fmin(self.ws,np.float64(fmin))
        
        if fatol is not None:
            optimpacklegacy.py_vmlmb_set_fatol(self.ws,np.float64(fatol))
        
        if frtol is not None:
            optimpacklegacy.py_vmlmb_set_frtol(self.ws,np.float64(frtol))
        
        if delta is not None:
            optimpacklegacy.py_vmlmb_set_delta(self.ws,np.float64(delta))
        
        if epsilon is not None:
            optimpacklegacy.py_vmlmb_set_epsilon(self.ws,np.float64(epsilon))
        
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
        return "Optimpack Optimizer [iter=%u] [task=%u] [status=%u]"%(self.iter,self.task,self.status)
    
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
        
        if self.verbose:
            print("OPL: %s"%optimpacklegacy.py_vmlmb_get_reason(self.ws))
        
        
    def run(self):
        if self.verbose:
            print("Optimizer starts")
        self.status = STATUS_NOT_FINISHED
        while (self.iter<self.itmax) and (self.task not in [TASK_CONV,TASK_ERROR,TASK_WARN]):
            self.step()
        
        message = optimpacklegacy.py_vmlmb_get_reason(self.ws)
        
        if self.iter>=self.itmax:
            self.status = STATUS_ITMAX
            message = "max number of iterations reached"
        
        if self.verbose:
            print("Optimizer ends (%s)"%message)
        
        return self.status

    
    def restart(self):
        optimpacklegacy.py_vmlmb_restart(self.ws)    
    
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
		
        if not type(self.x) is np.ndarray:
            raise TypeError("`x` must be a numpy.ndarray")

        if type(f) is np.float64:
            self.func = decorate_array(self.func)
            f = self.func(self.x,**self.kwargs)

        if not type(f) is np.ndarray:
            raise TypeError("`func` must return a numpy.float64 or a numpy.ndarray of length=1")
		
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
            
        self.last_f = f #save func
        self.last_g = g #save grad

