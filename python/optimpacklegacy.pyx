#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 08 17:18:00 2019

Cython file to connect C optimpack library to python

@author: rfetick
"""

import cython
import numpy as np
cimport numpy as np

#### ALL CDEF ####

cdef extern from "optimpacklegacy.h":
    int opl_vmlmb_monolithic_workspace_size(int n, int m)
	
cdef extern from "optimpacklegacy.h":
    int opl_vmlmb_monolithic_workspace_init(char* buffer, int n, int m)

cdef extern from "optimpacklegacy.h":
    int opl_vmlmb_iterate(int ws, double x[], double *f, double g[], int isfree[], const double h[])

cdef extern from "optimpacklegacy.h":
    int opl_vmlmb_set_fmin(int ws, double value) 

cdef extern from "optimpacklegacy.h":
    int opl_vmlmb_set_fatol(int ws, double value) 

cdef extern from "optimpacklegacy.h":
    int opl_vmlmb_set_frtol(int ws, double value) 

cdef extern from "optimpacklegacy.h":
    int opl_vmlmb_set_delta(int ws, double value) 

cdef extern from "optimpacklegacy.h":
    int opl_vmlmb_set_epsilon(int ws, double value) 

cdef extern from "optimpacklegacy.h":
    int opl_vmlmb_restart(int ws)

cdef extern from "optimpacklegacy.h":
    char* opl_vmlmb_get_reason(int ws)

#### ALL DEF ####

def py_vmlmb_workspace_size(n: int, m: int) -> int:
    return opl_vmlmb_monolithic_workspace_size(n, m)

@cython.boundscheck(False)
@cython.wraparound(False) 
def py_vmlmb_monolithic_workspace_init(np.ndarray[char, ndim=1, mode="c"] buffer not None, n: int, m: int) -> int:
    return opl_vmlmb_monolithic_workspace_init(&buffer[0], n, m)

@cython.boundscheck(False)
@cython.wraparound(False) 
def py_vmlmb_iterate(int ws, np.ndarray[double, ndim=1, mode="c"] x, np.ndarray[double, ndim=1, mode="c"] f, np.ndarray[double, ndim=1, mode="c"] g,np.ndarray[int, ndim=1, mode="c"] isfree)->int:
    return opl_vmlmb_iterate(ws, &x[0], &f[0], &g[0], &isfree[0], NULL)

def py_vmlmb_set_fmin(int ws, double value)->int:
    return opl_vmlmb_set_fmin(ws, value)

def py_vmlmb_set_fatol(int ws, double value)->int:
    return opl_vmlmb_set_fatol(ws, value)

def py_vmlmb_set_frtol(int ws, double value)->int:
    return opl_vmlmb_set_frtol(ws, value)

def py_vmlmb_set_delta(int ws, double value)->int:
    return opl_vmlmb_set_delta(ws, value)

def py_vmlmb_set_epsilon(int ws, double value)->int:
    return opl_vmlmb_set_epsilon(ws, value)

def py_vmlmb_restart(int ws)->int:
    return opl_vmlmb_restart(ws)
    
@cython.boundscheck(False)
@cython.wraparound(False)
def py_vmlmb_get_reason(int ws)->str:
    return opl_vmlmb_get_reason(ws).decode("utf-8")

