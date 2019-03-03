# -*- coding: utf-8; mode: python -*-
"""
Cython file defining a low-level Python interface to the OptimPacklegacy
library.

@author: rfetick, emmt
"""

cdef extern from "../src/optimpacklegacy.h":
    # Status returned in case of success, any other value is an error code.
    enum: OPL_SUCCESS
    enum: OPL_STP_EQ_STPMIN
    enum: OPL_STP_EQ_STPMAX
    enum: OPL_XTOL_TEST_SATISFIED
    enum: OPL_ROUNDING_ERROR
    enum: OPL_STPMAX_LT_STPMIN
    enum: OPL_STPMIN_LT_ZERO
    enum: OPL_XTOL_LT_ZERO
    enum: OPL_FTOL_LE_ZERO
    enum: OPL_GTOL_LE_ZERO
    enum: OPL_NOT_A_DESCENT
    enum: OPL_STP_GT_STPMAX
    enum: OPL_STP_LT_STPMIN
    enum: OPL_F_LE_FMIN
    enum: OPL_NOT_POSITIVE_DEFINITE
    enum: OPL_INSUFFICIENT_MEMORY
    enum: OPL_ILLEGAL_ADDRESS
    enum: OPL_INVALID_ARGUMENT
    enum: OPL_OUT_OF_BOUNDS
    enum: OPL_CORRUPTED
    enum: OPL_OVERFLOW
    enum: OPL_SYSTEM_ERROR

    # Possible values for an optimization task.
    enum: OPL_TASK_START    # start line search
    enum: OPL_TASK_FG       # caller has to compute function and gradient
    enum: OPL_TASK_FREEVARS # caller has to determine the free variables
    enum: OPL_TASK_NEWX     # new variables available for inspection
    enum: OPL_TASK_CONV     # search has converged
    enum: OPL_TASK_WARN     # search aborted with warning
    enum: OPL_TASK_ERROR    # search aborted with error

    # Opaque structure for VMLM-B workspace.
    ctypedef struct opl_vmlmb_workspace_t:
        pass

    # Prototype of called functions.
    const char* opl_get_default_message(int status)
    opl_vmlmb_workspace_t* opl_vmlmb_create(long n, long m)
    void opl_vmlmb_destroy(opl_vmlmb_workspace_t* ws)
    opl_vmlmb_workspace_t* opl_vmlmb_set_defaults(opl_vmlmb_workspace_t* ws)
    size_t opl_vmlmb_monolithic_workspace_size(long n, long m)
    opl_vmlmb_workspace_t* opl_vmlmb_monolithic_workspace_init(void* buf,
                                                               long n, long m)
    int opl_vmlmb_iterate(opl_vmlmb_workspace_t* ws, double x[], double *f,
                          double g[], int isfree[], const double h[])
    int opl_vmlmb_restart(opl_vmlmb_workspace_t* ws)
    int opl_vmlmb_restore(opl_vmlmb_workspace_t* ws,
                          double x[], double *f, double g[])
    int opl_vmlmb_set_fmin(opl_vmlmb_workspace_t* ws, double value)
    int opl_vmlmb_set_fatol(opl_vmlmb_workspace_t* ws, double value)
    int opl_vmlmb_set_frtol(opl_vmlmb_workspace_t* ws, double value)
    int opl_vmlmb_set_delta(opl_vmlmb_workspace_t* ws, double value)
    int opl_vmlmb_set_epsilon(opl_vmlmb_workspace_t* ws, double value)
    int opl_vmlmb_set_sxtol(opl_vmlmb_workspace_t* ws, double value)
    int opl_vmlmb_set_sftol(opl_vmlmb_workspace_t* ws, double value)
    int opl_vmlmb_set_sgtol(opl_vmlmb_workspace_t* ws, double value)
    int opl_vmlmb_get_task(opl_vmlmb_workspace_t* ws)
    int opl_vmlmb_get_status(opl_vmlmb_workspace_t* ws)
    const char* opl_vmlmb_get_reason(opl_vmlmb_workspace_t* ws)
    double opl_vmlmb_get_fmin(opl_vmlmb_workspace_t* ws)
    double opl_vmlmb_get_sftol(opl_vmlmb_workspace_t* ws)
    double opl_vmlmb_get_sgtol(opl_vmlmb_workspace_t* ws)
    double opl_vmlmb_get_sxtol(opl_vmlmb_workspace_t* ws)
    double opl_vmlmb_get_frtol(opl_vmlmb_workspace_t* ws)
    double opl_vmlmb_get_fatol(opl_vmlmb_workspace_t* ws)
    double opl_vmlmb_get_epsilon(opl_vmlmb_workspace_t* ws)
    double opl_vmlmb_get_delta(opl_vmlmb_workspace_t* ws)
    double opl_vmlmb_get_step(opl_vmlmb_workspace_t* ws)
    double opl_vmlmb_get_gnorm(opl_vmlmb_workspace_t* ws)
    long opl_vmlmb_get_n(opl_vmlmb_workspace_t* ws)
    long opl_vmlmb_get_m(opl_vmlmb_workspace_t* ws)
    long opl_vmlmb_get_evaluations(opl_vmlmb_workspace_t* ws)
    long opl_vmlmb_get_iterations(opl_vmlmb_workspace_t* ws)
    long opl_vmlmb_get_restarts(opl_vmlmb_workspace_t* ws)
