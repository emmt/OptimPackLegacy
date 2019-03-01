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
    """`clamp(x, lo = None, hi = None)` yields the result of clamping the values
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

def zipfg(func, grad):
    """
    `zipfg(func,grad)` yields a function which, when called with argument `x`,
    returns the tuple `(func(x),grad(x))`.  This may be useful to combine an
    objective function and its gradient for VMLM-B algorithm.  However note
    that it is usually more efficient to compute an objective function and its
    gradient in a same function because they involve similar terms."""
    return lambda x: (func(x), grad(x))

try:
    import resource
    def cputime():
        """cputime()

        yields the user CPU time in seconds used since the start of the
        process."""
        return resource.getrusage(resource.RUSAGE_SELF)[0]

except ImportError:
    import time
    # There is no distinction of user/system time under Windows, so we just use
    # time.clock() with respect to some origin...
    def cputime():
        """cputime()

        yields the time in seconds since the Epoch.  This is because, under
        Windows, the consumed CPU time can't be measured."""
        return time.clock()

# Default printer for VMLM-B algorithm.
def _default_vmlmb_printer(opt, x, f, cpu):
    neval = opt.get_evaluations()
    if neval == 0:
        print("ITER  EVAL   CPU (ms)        FUNC               GNORM   STEPLEN")
        print("---------------------------------------------------------------")
    else:
        niter = opt.get_iterations()
        gnorm = opt.get_gnorm()
        stpln = opt.get_step()
        print("%4d %5d %10.3f  %+-24.15e%-9.1e%-9.1e" % (niter, neval,
                                                         cpu*1e-3, f,
                                                         gnorm, stpln))

def vmlmb(fg, x, xmin=None, xmax=None, mem=None, h=None, fmin=None,
          verb=-1, quiet=False, printer=_default_vmlmb_printer,
          maxiter=None, maxeval=None, output=None,
          frtol=None, fatol=None, gatol=None, grtol=None,
          sftol=None, sgtol=None, sxtol=None):
    """vmlmb(fg, x0, xmin=None, xmax=None, mem=None, h=None, fmin=None,
          verb=-1, quiet=False, viewer=None, printer=None,
          maxiter=None, maxeval=None, output=None,
          frtol=None, fatol=None, gatol=None, grtol=None,
          sftol=None, sgtol=None, sxtol=None) -> (x, fx, gx)

    Find a minimum of a multivariate function by an iterative minimization
    algorithm (a limited memory quasi-Newton method) possibly with simple bound
    constraints on the parameters.


    Result
    ------

    The solution `x`, the corresponding function value `fx` and gradient `gx`
    are returned.  In case of early termination, the best solution found so far
    is returned.  If the multivariate function has more than one minimum, which
    minimum is returned is undefined (although it depends on the starting
    variables `x0`).


    Parameters
    ----------

    `fg` - User defined function implementing the computation of the objective
         function value and its gradient.  Its prototype is:

         def fg(x):
              fx = .... # compute function value at X
              gx = .... # store gradient of F in GX
              return (fx, gx)

         This is motivated by the fact that it is usually more efficient to
         compute an objective function and its gradient in a same function
         because they involve similar terms.  If the objective function and its
         gradient are however implemented by two different functions, say
         `func` and `grad`, `zipfg(func,grad)` can be used to specify `fg`.

     `x0` - The starting solution (an array).

     `xmin`, `xmax` - Lower and upper bounds for `x`.  Must be conformable with
         the sought variables `x`.  For instance with `xmin=0`, the
         non-negative solution will be returned.  Not taken into account when
         set to `None` (the default).

     `mem` - Number of previous iterates to memorize to compute the LBFGS
         approximation of the inverse Hessian matrix (default:
         `mem=min(5,len(x0))`).

     `h` - Diagonal preconditioner.

     `maxiter` - Maximum number of iterations (default: no limits).

     `maxeval` - Maximum number of function evaluations (default: no limits).

     `verb` - Must be an integer.  If greater or equal one, print out
         information every `verb` iterations and for the final one.

     `quiet` - If true and not in verbose mode, do not print warning nor
         convergence error messages.

     `printer` - User defined function to call every `verb` iterations (see
         argument `verb` above) to printout iteration information.  The
         function will be called as:

            printer(opt, x, f, cpu)

         where `opt` is the VMLM-B workspace used to perform the iterations,
         `x` is the current solution, `f` is the function value and `cpu` is
         the elapsed CPU time in seconds.  Methods `opt.get_evaluations()` or
         `opt.get_iterations()` can be used to retrieve the number of
         evaluations of the objective function and the number of iterations of
         the algorithm.  If the number of evaluations is equal to zero, the
         algorithm has not yet started and only some head lines shall be
         printed (and nothing else).  Otherwise, information such as the
         Euclidean norm of the gradient and the length of the step along the
         search direction can be retrieved with `opt.get_gnorm()` or
         `opt.get_step()`.

     `fatol`, `frtol` - Relative function change tolerance for convergence
         (default: 1.5e-8).

     `gatol`, `grtol` - Absolute and relative gradient tolerances for
         convergence which is assumed whenever the Euclidean norm of the
         (projected) gradient is smaller or equal `max(gatol,grtol*ginit)`
         where `ginit` is the Euclidean norm of the (projected) initial
         gradient.  By default, `gatol=0` and `grtol=1e-6`.

     `sftol`, `sgtol`, `sxtol` - Line search tolerance and safeguard
         parameters.

    """

    # Starting variables.
    x = np.array(x, copy=True, dtype=np.double, order='C')
    f = np.empty(1, dtype=np.double, order='C')
    #g = np.empty(x.shape, dtype=np.double, order='C')
    freevars = None

    # Create optimizer workspace with options.
    opt = VMLMB(len(x), mem)
    if fmin is not None: opt.set_fmin(fmin)
    if frtol is not None: opt.set_frtol(frtol)
    if fatol is not None: opt.set_fatol(fatol)
    if grtol is not None: opt.set_grtol(grtol)
    if gatol is not None: opt.set_gatol(gatol)
    if sftol is not None: opt.set_sftol(sftol)
    if sgtol is not None: opt.set_sgtol(sgtol)
    if sxtol is not None: opt.set_sxtol(sxtol)

    # Start iterations.
    neval = 0 # number of objective function evaluations -- also given by
              # opt.get_evaluations()
    niter = 0 # number of algorithm iterations -- also given by
              # opt.get_iterations()
    task = opt.get_task()
    stop = False
    mesg = None
    if verb > 0:
        t0 = cputime()
        printer(opt, x, 0.0, 0.0)
    while True:
        prnt = False
        if task == TASK_FG:
            if maxeval is not None and neval >= maxeval:
                # Too many objective function evaluations.  We restore the
                # variables at the start of the line search which is a cheap
                # way (no extra memory cost) to recover variables which should
                # be nearly the best ones.  We wait until last iteration
                # information have been printed to stop the iterations.
                mesg = 'too many function evaluations ({})'.format(neval)
                stop = True
                prnt = (verb > 0)
                opt.restore(x, f, g)
            else:
                # Apply the bounds (if any), then compute the objective
                # function and its gradient.
                x = clamp(x, xmin, xmax)
                (f[0],g) = fg(x)
                neval += 1
        elif task == TASK_FREEVARS:
            # Determine the set of unbind variables.
            freevars = get_free_variables(x, g, xmin, xmax)
        elif task == TASK_NEWX:
            # Line search has converged or objective function has been computed
            # for the initial variables.  We may be interested in printing some
            # information about the current iterate.
            niter = opt.get_iterations()
            if maxiter is not None and niter >= maxiter:
                # Too many algorithm iterations.  See comments above in the
                # case of too many evaluations.
                mesg = 'too many iterations ({})'.format(niter)
                stop = True
                prnt = (verb > 0)
            else:
                prnt = (verb > 0 and (niter % verb) == 0)
        else:
            if task == TASK_CONV:
                mesg = 'algorithm has converged'
            elif task == TASK_WARN:
                mesg = 'algorithm has warnings: ' + opt.get_reason()
                opt.restore(x, f, g) # FIXME: is this needed?
            elif task == TASK_ERROR:
                # In case of more serious errors, we do not restore the best
                # solution so far because it may be useful to examine the
                # current variables to figure out the origin of the problem.
                mesg = 'algorithm has errors: ' + opt.get_reason()
            else:
                mesg = 'illegal task value: ' + repr(task)
            stop = True
            prnt = (verb > 0)

        if prnt:
            printer(opt, x, f[0], cputime() - t0)
        if stop:
            if not quiet: print(mesg)
            break

        # Perform next algorithm step.
        task = opt.iterate(x, f, g, freevars, h)

    # Return the solution so far, the function value and the gradient.
    return (x, f[0], g)

#------------------------------------------------------------------------------
# The following class encapsulates the low level interface to the
#  VMLM-B algorithm in OptimPackLegacy library.

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
