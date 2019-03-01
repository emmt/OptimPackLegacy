# -*- coding: utf-8; mode: python -*-
"""
Cython file to provide a Python interface to the OptimPackLegacy library.

@author: rfetick, emmt
"""

import cython
import numpy as np
cimport numpy as np

# Import low-level interface to C library (defined in coptimpacklegacy.pxd).
cimport coptimpacklegacy as copl

# NumPy type for logical values in OptimPackLegacy.  FIXME: This type should be
# automatically determined.
LOGICAL = np.int32

# Make status codes available.  Any other value thena SUCCESS is an error code.
SUCCESS               = copl.OPL_SUCCESS
STP_EQ_STPMIN         = copl.OPL_STP_EQ_STPMIN
STP_EQ_STPMAX         = copl.OPL_STP_EQ_STPMAX
XTOL_TEST_SATISFIED   = copl.OPL_XTOL_TEST_SATISFIED
ROUNDING_ERROR        = copl.OPL_ROUNDING_ERROR
STPMAX_LT_STPMIN      = copl.OPL_STPMAX_LT_STPMIN
STPMIN_LT_ZERO        = copl.OPL_STPMIN_LT_ZERO
XTOL_LT_ZERO          = copl.OPL_XTOL_LT_ZERO
FTOL_LE_ZERO          = copl.OPL_FTOL_LE_ZERO
GTOL_LE_ZERO          = copl.OPL_GTOL_LE_ZERO
NOT_A_DESCENT         = copl.OPL_NOT_A_DESCENT
STP_GT_STPMAX         = copl.OPL_STP_GT_STPMAX
STP_LT_STPMIN         = copl.OPL_STP_LT_STPMIN
F_LE_FMIN             = copl.OPL_F_LE_FMIN
NOT_POSITIVE_DEFINITE = copl.OPL_NOT_POSITIVE_DEFINITE
INSUFFICIENT_MEMORY   = copl.OPL_INSUFFICIENT_MEMORY
ILLEGAL_ADDRESS       = copl.OPL_ILLEGAL_ADDRESS
INVALID_ARGUMENT      = copl.OPL_INVALID_ARGUMENT
OUT_OF_BOUNDS         = copl.OPL_OUT_OF_BOUNDS
CORRUPTED             = copl.OPL_CORRUPTED
OVERFLOW              = copl.OPL_OVERFLOW
SYSTEM_ERROR          = copl.OPL_SYSTEM_ERROR

# Make task codes available for others.
TASK_START    = copl.OPL_TASK_START    # start line search
TASK_FG       = copl.OPL_TASK_FG       # caller has to compute function and gradient
TASK_FREEVARS = copl.OPL_TASK_FREEVARS # caller has to determine the free variables
TASK_NEWX     = copl.OPL_TASK_NEWX     # new variables available for inspection
TASK_CONV     = copl.OPL_TASK_CONV     # search has converged
TASK_WARN     = copl.OPL_TASK_WARN     # search aborted with warning
TASK_ERROR    = copl.OPL_TASK_ERROR    # search aborted with error

class Error(Exception):
    """Raised when an error occurs in a call to a fcuntion of the OptimPackLegacy
    library.  The associated value is an integer status code."""
    def __init__(self, code):
        self.code = code
    def __str__(self):
        return get_message(self.code) + " [" + repr(self.code) + "]"

class SizeError(Exception):
    """Raised when an argument has the bad number of elements.  The associated
    values are a textual description of the argument and the expected number of
    elements."""
    def __init__(self, attr, size):
        self.attr = attr
        self.size = size
    def __str__(self):
        return self.attr + " must be a NumPy array of length " + repr(self.size)

def get_message(code) -> str:
    """`get_message(code)` yields the textual message corresponding to the
    value of `code` (an integer returned by one of the functions of the
    OptimPackLegacy library."""
    return copl.opl_get_default_message(code)#.decode("utf-8")

def clamp(x, lo = None, hi = None):
    """
    `clamp(x, lo = None, hi = None)` yields the result of clamping the values
    of array `x` between the lower and upper bounds respectively specified by
    `lo` and `hi`.  A given bound is ignored if the corresponding argument is
    `None`."""
    if lo is not None:
        x = np.maximum(x, lo)
    if hi is None:
        return x
    else:
        return np.minimum(x, hi)

def get_free_variables(x, g, lo = None, hi = None, dest = None):
    if dest is None:
        dest = np.ones(x.shape, dtype=LOGICAL, order='C')
        #dest = np.empty(x.shape, dtype=LOGICAL, order='C')
        #dest[:] = 1
    else:
        dest[:] = 1
    if lo is not None:
        # The following expression is the same as performing the elementwise
        # operation x > lo || g < 0.
        dest *= np.maximum(x > lo, g < 0).astype(LOGICAL)
    if hi is not None:
        # The following expression is the same as performing the elementwise
        # operation x < hi || g > 0.
        dest *= np.maximum(x < hi, g > 0).astype(LOGICAL)
    return dest

cdef class VMLMB:
    """
    `VMLMB(n, m)` yields a new workspace for solving an optimization problem in
    `n` variables by means of a limited-memory quasi-Newton method with `m`
    memorized previous iterates.  Simple bound constraints may be applied to
    the variables.  Algorithm is VMLM-B which expects the objective function
    and its gradient.
    """

    # Pointer to the opaque workspace structure.
    cdef copl.opl_vmlmb_workspace_t* _ws


    # Allocate a new VMLM-B workspace when object is created.
    def __cinit__(self, n, mem = None):
        if mem is None:
            mem = min(n, 5)
        self._ws = copl.opl_vmlmb_create(n, mem)
        if self._ws is NULL:
            raise MemoryError()

    # Manage to free allocated workspace when object is finalized.
    def __dealloc__(self):
        if self._ws is not NULL:
            copl.opl_vmlmb_destroy(self._ws)

    def restart(self):
        status = copl.opl_vmlmb_restart(self._ws)
        if status != SUCCESS:
            raise Error(status)

    def get_fmin(self):
        return copl.opl_vmlmb_get_fmin(self._ws)

    def set_fmin(self, value):
        status = copl.opl_vmlmb_set_fmin(self._ws, value)
        if status != SUCCESS:
            raise Error(status)

    def get_sftol(self):
        return copl.opl_vmlmb_get_sftol(self._ws)

    def set_sftol(self, value):
        status = copl.opl_vmlmb_set_sftol(self._ws, value)
        if status != SUCCESS:
            raise Error(status)

    def get_sgtol(self):
        return copl.opl_vmlmb_get_sgtol(self._ws)

    def set_sgtol(self, value):
        status = copl.opl_vmlmb_set_sgtol(self._ws, value)
        if status != SUCCESS:
            raise Error(status)

    def get_sxtol(self):
        return copl.opl_vmlmb_get_sxtol(self._ws)

    def set_sxtol(self, value):
        status = copl.opl_vmlmb_set_sxtol(self._ws, value)
        if status != SUCCESS:
            raise Error(status)

    def get_frtol(self):
        return copl.opl_vmlmb_get_frtol(self._ws)

    def set_frtol(self, value):
        status = copl.opl_vmlmb_set_frtol(self._ws, value)
        if status != SUCCESS:
            raise Error(status)

    def get_fatol(self):
        return copl.opl_vmlmb_get_fatol(self._ws)

    def set_fatol(self, value):
        status = copl.opl_vmlmb_set_fatol(self._ws, value)
        if status != SUCCESS:
            raise Error(status)

    def get_epsilon(self):
        return copl.opl_vmlmb_get_epsilon(self._ws)

    def set_epsilon(self, value):
        status = copl.opl_vmlmb_set_epsilon(self._ws, value)
        if status != SUCCESS:
            raise Error(status)

    def get_delta(self):
        return copl.opl_vmlmb_get_delta(self._ws)

    def set_delta(self, value):
        status = copl.opl_vmlmb_set_delta(self._ws, value)
        if status != SUCCESS:
            raise Error(status)

    def get_status(self):
        """Yield the current status."""
        return copl.opl_vmlmb_get_status(self._ws)

    def get_reason(self) -> str:
        """Yield a textual explanation of the next task to perform or of the
        current error if any."""
        return copl.opl_vmlmb_get_reason(self._ws)#.decode("utf-8")

    def get_task(self):
        """Yield the next task to perform."""
        return copl.opl_vmlmb_get_task(self._ws)

    def get_step(self):
        """Yield the length of the step along the search direction."""
        return copl.opl_vmlmb_get_step(self._ws)

    def get_gnorm(self):
        """Yield the Euclidean norm of the gradient."""
        return copl.opl_vmlmb_get_gnorm(self._ws)

    def get_size(self):
        """Yield the number of variables."""
        return copl.opl_vmlmb_get_n(self._ws)

    def get_memory(self):
        """Yield the maximum number of memorized previous iterates."""
        return copl.opl_vmlmb_get_m(self._ws)

    def get_evaluations(self):
        """Yield the number of objective function evaluations."""
        return copl.opl_vmlmb_get_evaluations(self._ws)

    def get_iterations(self):
        """Yield the number of algorithm iterations."""
        return copl.opl_vmlmb_get_iterations(self._ws)

    def get_restarts(self):
        """Yield the number of algorithm restarts."""
        return copl.opl_vmlmb_get_restarts(self._ws)

    # Turn off bounds checking and wrapping of negative indices.
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def iterate(self,
                np.ndarray[double, ndim=1, mode='c'] x,
                np.ndarray[double, ndim=1, mode='c'] f,
                np.ndarray[double, ndim=1, mode='c'] g,
                np.ndarray[int,    ndim=1, mode='c'] isfree = None,
                np.ndarray[double, ndim=1, mode='c'] h = None):
        # Check sizes of arguments.
        n = self.get_size()
        if len(f) != 1: raise SizeError("objective function value `f`", 1)
        if len(x) != n: raise SizeError("variables `x`", n)
        if len(g) != n: raise SizeError("gradient `g`", n)
        if isfree is not None and len(isfree) != n:
            raise SizeError("mask of free variables `isfree`", n)
        if h is not None and len(h) != n:
            raise SizeError("diagonal preconditioner `h`", n)

        # Mask of free variables and diagonal preconditioner are optional.
        if isfree is None:
            # Mask of free variables has not been specified (or is None)
            # meaning that there are no bound constraints.
            if h is None:
                task = copl.opl_vmlmb_iterate(self._ws, &x[0], &f[0], &g[0],
                                              NULL, NULL)
            else:
                task = copl.opl_vmlmb_iterate(self._ws, &x[0], &f[0], &g[0],
                                              NULL, &h[0])
        else:
            # Mask of free variables has been specified to account for bound
            # constraints.
            if h is None:
                task = copl.opl_vmlmb_iterate(self._ws, &x[0], &f[0], &g[0],
                                              &isfree[0], NULL)
            else:
                task = copl.opl_vmlmb_iterate(self._ws, &x[0], &f[0], &g[0],
                                              &isfree[0], &h[0])
        return task

    # Turn off bounds checking and wrapping of negative indices.
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def restore(self,
                np.ndarray[double, ndim=1, mode='c'] x,
                np.ndarray[double, ndim=1, mode='c'] f,
                np.ndarray[double, ndim=1, mode='c'] g):
        # Check sizes of arguments and call low-level function.
        n = self.get_size()
        if len(f) != 1: raise SizeError("objective function value `f`", 1)
        if len(x) != n: raise SizeError("variables `x`", n)
        if len(g) != n: raise SizeError("gradient `g`", n)
        return copl.opl_vmlmb_restore(self._ws, &x[0], &f[0], &g[0])
