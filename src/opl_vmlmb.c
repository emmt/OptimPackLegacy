/*
 * opl_vmlmb.c --
 *
 * Variable Metric Limited Memory with Bound constraints for OptimPackLegacy
 * library.
 *
 *-----------------------------------------------------------------------------
 *
 * This file is part of OptimPackLegacy
 * <https://github.com/emmt/OptimPackLegacy>.
 *
 * Copyright (c) 2003-2021, Éric Thiébaut.
 *
 * OptimPackLegacy is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option)
 * any later version.
 *
 * OptimPackLegacy is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * OptimPackLegacy (file "LICENSE" in the top source directory); if not, write
 * to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA 02111-1307 USA
 *
 *-----------------------------------------------------------------------------
 */

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>
#include "opl_private.h"

typedef unsigned char byte_t;

/*
 * Some constants.
 *
 * See http://stackoverflow.com/questions/1923837/how-to-use-nan-and-inf-in-c
 * for a discussion about how to define NaN and infinite.  Alternatives for Inf
 * and NaN:
 *
 *   double NaN = strtod("NaN", NULL);
 *   double Inf = strtod("Inf", NULL);
 *
 * but cannot be constants then.  Some macros may be defined in <math.h> (as
 * part of C99 standard):
 *
 *   #include <math.h>
 *   #ifdef NAN
 *   const double NaN = NAN;
 *   #end
 *   #ifdef INFINITY
 *   const double Inf = INFINITY;
 *   #end
 */
#if __STDC_VERSION__ >= 199901L /* C99 */
#  define CHOICE_(a, b) a
#else                       /* not C99 */
#  define CHOICE_(a, b) b
#endif
static const double NaN = CHOICE_(NAN, 0.0/0.0);
#if 0 /* not used */
static const double Inf = CHOICE_(INFINITY, 1.0/0.0);
#endif

/* Default settings. */
#define DEFAULT_STPMAX   1e20
#define DEFAULT_SFTOL    0.001
#define DEFAULT_SGTOL    0.9
#define DEFAULT_SXTOL    0.1
#define DEFAULT_FATOL    0.0
#define DEFAULT_FRTOL    1e-10
#define DEFAULT_DELTA    1e-3
#define DEFAULT_EPSILON  0.0

#define FLAG_FMIN        (1 << 0)

/* `task`, `context` and some other members are shared with the embedded
   line search structure. */
#define TASK(ws)    (ws)->lnsrch.task
#define CONTEXT(ws) (ws)->lnsrch.context
#define STATUS(ws)  (ws)->lnsrch.context.status
#define SXTOL(ws)   (ws)->lnsrch.xtol
#define SFTOL(ws)   (ws)->lnsrch.ftol
#define SGTOL(ws)   (ws)->lnsrch.gtol
#define SXTOL(ws)   (ws)->lnsrch.xtol
#define STPMIN(ws)  (ws)->lnsrch.stpmin
#define STPMAX(ws)  (ws)->lnsrch.stpmax

/* Allocate an array of given type and size. */
#define NEW(type, number) ((type*)malloc((number)*sizeof(type)))

/*---------------------------------------------------------------------------*/
/* PRIVATE FUNCTIONS */

static double min(double a, double b) {
    return (a <= b ? a : b);
}

static double max(double a, double b) {
    return (a >= b ? a : b);
}

/* Apply L-BFGS recursion to compute search direction. */
static void compute_direction(
    opl_vmlmb_workspace_t* ws,
    double d[],
    const int isfree[],
    const double h[]);

/* Compute next step and set task. */
static opl_task_t next_step(
    opl_vmlmb_workspace_t* ws,
    double x[]);

/* Check H and ISFREE arrays, possibly fix ISFREE, returns non-zero
   in case of error. */
static int check_free(
    opl_vmlmb_workspace_t* ws,
    int isfree[],
    const double h[]);

/* Report a failure. */
static opl_task_t failure(
    opl_vmlmb_workspace_t* ws,
    opl_status_t status,
    const char* mesg)
{
    opl_set_context(&CONTEXT(ws), status, mesg, OPL_PERMANENT);
    return (TASK(ws) = OPL_TASK_ERROR);
}

/* Report a success. */
static opl_task_t success(
    opl_vmlmb_workspace_t* ws,
    opl_task_t task,
    const char* mesg)
{
    opl_set_context(&CONTEXT(ws), OPL_SUCCESS, mesg, OPL_PERMANENT);
    return (TASK(ws) = task);
}

/*---------------------------------------------------------------------------*/

#define GET_MEMBER(type, name, rvalue, invalid) \
    type opl_vmlmb_get_##name(                  \
        opl_vmlmb_workspace_t* ws)              \
    {                                           \
        if (ws == NULL) {                       \
            errno = EFAULT;                     \
            return invalid;                     \
        }                                       \
        return rvalue;                          \
    }

GET_MEMBER(double, fmin, ws->fmin, NaN);
GET_MEMBER(double, fatol, ws->fatol, NaN);
GET_MEMBER(double, frtol, ws->frtol, NaN);
GET_MEMBER(double, delta, ws->delta, NaN);
GET_MEMBER(double, epsilon, ws->epsilon, NaN);
GET_MEMBER(double, sxtol, SXTOL(ws), NaN);
GET_MEMBER(double, sftol, SFTOL(ws), NaN);
GET_MEMBER(double, sgtol, SGTOL(ws), NaN);
GET_MEMBER(double, step, ws->stp, NaN);
GET_MEMBER(double, gnorm, ws->gnorm, NaN);
GET_MEMBER(long, n, ws->n, 0);
GET_MEMBER(long, m, ws->m, 0);
GET_MEMBER(opl_task_t, task, TASK(ws), OPL_TASK_ERROR);
GET_MEMBER(opl_status_t, status, STATUS(ws), OPL_ILLEGAL_ADDRESS);
GET_MEMBER(long, iterations, ws->iterations, -1);
GET_MEMBER(long, evaluations, ws->evaluations, -1);
GET_MEMBER(long, restarts, ws->restarts, -1);
GET_MEMBER(const char*, reason, opl_get_message(&CONTEXT(ws)), NULL);

#undef GET_MEMBER

opl_status_t opl_vmlmb_set_fmin(
    opl_vmlmb_workspace_t* ws,
    double value)
{
    if (ws == NULL) {
        errno = EFAULT;
        return OPL_ILLEGAL_ADDRESS;
    }
    if (isnan(value) || value < -DBL_MAX) {
        ws->flags &= ~FLAG_FMIN;
        ws->fmin = NaN;
    } else {
        ws->fmin = value;
    }
    return OPL_SUCCESS;
}

#define SET_MEMBER(name, lvalue, invalid)       \
    opl_status_t opl_vmlmb_set_##name(          \
        opl_vmlmb_workspace_t* ws,              \
        double value)                           \
    {                                           \
        if (ws == NULL) {                       \
            errno = EFAULT;                     \
            return OPL_ILLEGAL_ADDRESS;         \
        }                                       \
        if (invalid) {                          \
            errno = EINVAL;                     \
            return OPL_INVALID_ARGUMENT;        \
        }                                       \
        lvalue = value;                         \
        return OPL_SUCCESS;                     \
    }

SET_MEMBER(fatol, ws->fatol, value < 0);
SET_MEMBER(frtol, ws->frtol, value < 0);
SET_MEMBER(delta, ws->delta, value < 0);
SET_MEMBER(epsilon, ws->epsilon, value < 0);
SET_MEMBER(sxtol, SXTOL(ws), value <= 0 || value >= 1);
SET_MEMBER(sftol, SFTOL(ws), value <= 0 || value >= 1);
SET_MEMBER(sgtol, SGTOL(ws), value <= 0 || value >= 1);

#undef SET_MEMBER

opl_task_t opl_vmlmb_restart(
    opl_vmlmb_workspace_t* ws)
{
    if (ws == NULL) {
        errno = EFAULT;
        return OPL_TASK_ERROR;
    }
    ws->evaluations = 0;
    ws->iterations = 0;
    ws->restarts = 0;
    ws->mark = -1;
    ws->mp = 0;
    ws->searching = OPL_FALSE;
    ws->f0 = 0;
    ws->gd = 0;
    ws->g0d = 0;
    ws->stp = 0;
    STPMIN(ws) = 0;
    STPMAX(ws) = 0;
    ws->gnorm = -1;
    ws->g0norm = -1;
    return success(ws, OPL_TASK_FG, "compute f(x) and g(x)");
}

opl_task_t opl_vmlmb_warm_restart(
    opl_vmlmb_workspace_t* ws)
{
    if (ws == NULL) {
        errno = EFAULT;
        return OPL_TASK_ERROR;
    }
    ws->searching = OPL_FALSE;
    return success(ws, OPL_TASK_FG, "compute f(x) and g(x)");
}

opl_task_t opl_vmlmb_restore(
    opl_vmlmb_workspace_t* ws,
    double x[],
    double *f,
    double g[])
{
    if (ws == NULL || x == NULL || f == NULL || g == NULL) {
        errno = EFAULT;
        return OPL_TASK_ERROR;
    }
    if (TASK(ws) == OPL_TASK_FG && ws->evaluations > 0) {
        *f = ws->f0;
        ws->gnorm = ws->g0norm;
        opl_dcopy(ws->n, ws->S[ws->mark], x);
        opl_dcopy(ws->n, ws->Y[ws->mark], g);
        success(ws, OPL_TASK_NEWX,
                "restored solution available for inspection");
    }
    return TASK(ws);
}

/*---------------------------------------------------------------------------*/
/* PERFORM NEXT VMLMB STEP */

opl_task_t opl_vmlmb_iterate(
    opl_vmlmb_workspace_t* ws,
    double x[],
    double *f,
    double g[],
    int isfree[],
    const double h[])
{
    /* Local variables. */
    double stpmax;
    long m = ws->m;
    long n = ws->n;
    double* d = ws->d;
    double** S_ = ws->S;
    double** Y_ = ws->Y;
    int have_fmin = ((ws->flags & FLAG_FMIN) != 0);

# define S(j) S_[j]
# define Y(j) Y_[j]

    switch (TASK(ws)) {

    case OPL_TASK_FG:
        /* Caller has performed a new evaluation of the function and its
           gradient. */
        ++ws->evaluations;
        if (have_fmin && *f <= ws->fmin) {
            return failure(ws, OPL_F_LE_FMIN, "initial F <= FMIN");
        }
        if (ws->searching) {
            /* A line seach is in progress.  Check for convergence of this
             * search. */
            opl_status_t status;
            opl_task_t task;

            /* Compute directional derivative for the line search. */
            ws->gd = -opl_ddot(n, g, d);

            /* Call line search iterator to check for line search
             * convergence. */
            opl_csrch_iterate(&ws->lnsrch, *f, ws->gd, &ws->stp);
            task = TASK(ws);
            if (task == OPL_TASK_FG) {
                /* Line search has not converged.  Compute a new iterate. */
                return next_step(ws, x);
            }
            status = STATUS(ws);
            if (task == OPL_TASK_CONV ||
                (task == OPL_TASK_WARN && (status == OPL_XTOL_TEST_SATISFIED ||
                                           status == OPL_ROUNDING_ERROR))) {
                /* Line search has converged. */
                ++ws->iterations;
                ws->searching = OPL_FALSE;
            } else {
                /* Error or warning in line search (task and context should
                   have been set so as to describe the problem).  Restore
                   solution at start of line search. */
                opl_dcopy(n, S(ws->mark), x);
                opl_dcopy(n, Y(ws->mark), g);
                ws->gnorm = ws->g0norm;
                *f = ws->f0;
                break;
            }
        }

        /* Request the caller to compute the set of unbinded variables. */
        return success(ws, OPL_TASK_FREEVARS, "determine free variables");

    case OPL_TASK_FREEVARS:
        /* Copy the (projected) gradient into D and compute the norm of the
           (projected) gradient. */
        if (check_free(ws, isfree, h) != 0) {
            break;
        }
        opl_dcopy_free(n, g, d, isfree);
        ws->gnorm = opl_dnrm2(n, d);
        if (ws->gnorm == 0) {
            return success(ws, OPL_TASK_CONV, "local minimum found");
        }
        if (ws->evaluations > 1) {
            /* Update L-BFGS model. (We always compute the effective step to
               account for bound constraints and, at least, numerical rounding
               or truncation errors.) */
            double* Smark = S(ws->mark);
            double sts = 0;
            for (long i = 0; i < n; ++i) {
                double s = Smark[i] - x[i];
                Smark[i] = s;
                sts += s*s;
            }
            double* Ymark = Y(ws->mark);
            double yty = 0;
            for (long i = 0; i < n; ++i) {
                double y = Ymark[i] - g[i];
                Ymark[i] = y;
                yty += y*y;
            }
            if (ws->mp < m) {
                ++ws->mp;
            }
            if (sts <= 0) {
                return success(ws, OPL_TASK_WARN, "no parameter change");
            }
            if (yty <= 0) {
                return success(ws, OPL_TASK_WARN, "no gradient change");
            }
            if (isfree == NULL) {
                double sty = opl_ddot(n, Smark, Ymark);
                ws->gamma = sty/yty;
                ws->rho[ws->mark] = sty;
            }

            /* Test for global convergence. */
            double fchange = max(fabs(*f - ws->f0), fabs(ws->stp*ws->g0d));
            if (fchange <= ws->frtol*fabs(ws->f0)) {
                return success(ws, OPL_TASK_CONV, "FRTOL test satisfied");
            }
            if (fchange <= ws->fatol) {
                return success(ws, OPL_TASK_CONV, "FATOL test satisfied");
            }
        }

        /* Set task to signal a new iterate. */
        return success(ws, OPL_TASK_NEWX,
                       "new improved solution available for inspection");

    case OPL_TASK_NEWX:
    case OPL_TASK_CONV:
    case OPL_TASK_WARN:
        /* Compute a new search direction.  D already contains the (projected)
           gradient. */
        if (ws->mp > 0) {
            /* Apply L-BFGS recursion to compute a search direction. */
            compute_direction(ws, d, isfree, h);
            if (ws->mp > 0) {
                /* Set initial step size and compute dot product of gradient
                   and search direction.  If L-BFGS recursion failed to produce
                   a sufficient descent search direction, the algorithm is
                   restarted with the steepest descent. */
                ws->stp = 1;
                ws->gd = -opl_ddot(n, g, d);
                if (ws->epsilon > 0 ?
                    (ws->gd > -ws->epsilon*ws->gnorm*opl_dnrm2(n, d)) :
                    (ws->gd >= 0)) {
                    /* Insufficient descent direction.  Manage to restart the
                       L-BFGS with steepest descent. */
                    ws->mp = 0;
                    opl_dcopy_free(n, g, d, isfree);
                }
            }
            if (ws->mp <= 0) {
                /* L-BFGS recursion has been restarted because it has failed to
                   produce an acceptable search direction. */
                ++ws->restarts;
            }
        }

        if (ws->mp == 0) {
            /* Compute initial search direction (or after a restart).  D is
               already the (projected) gradient. */
            if (h != NULL) {
                /* Use diagonal preconditioner to compute initial search
                   direction. */
                for (long i = 0; i < n; ++i) {
                    d[i] *= h[i];
                }
                ws->stp = 1;
                ws->gd = -opl_ddot(n, g, d);
                if (ws->gd >= 0) {
                    return failure(ws, OPL_NOT_POSITIVE_DEFINITE,
                                   "preconditioner is not positive definite");
                }
            } else {
                /* No preconditioning, use a small step along the steepest
                 * descent. */
                if (ws->delta > 0) {
                    ws->stp = (opl_dnrm2(n, x)/ws->gnorm)*ws->delta;
                } else {
                    ws->stp = 0;
                }
                if (ws->stp <= 0) {
                    /* The following step length is quite arbitrary but to
                       avoid this would require to know the typical size of the
                       variables. */
                    ws->stp = 1.0/ws->gnorm;
                }
                ws->gd = -ws->gnorm*ws->gnorm;
            }
        }

        /* Advance the mark and save point at start of line search.  Note that
           this point has forcibly been projected so it is feasible. */
        ws->mark = (ws->mark + 1)%m;
        ws->f0 = *f;
        ws->g0d = ws->gd;
        ws->g0norm = ws->gnorm;
        opl_dcopy(n, x, S(ws->mark)); /* save parameters X0 */
        opl_dcopy(n, g, Y(ws->mark)); /* save gradient G0 */

        /* Set step bounds and task to start line search. */
#if 1
        stpmax = DEFAULT_STPMAX;
#else
        if (have_fmin) {
            stpmax = (ws->fmin - ws->f0)/(SGTOL(ws)*ws->g0d);
        } else {
            double temp = fabs(ws->f0);
            if (temp < 1) {
                temp = 1;
            }
            stpmax = temp/(SGTOL(ws)*ws->g0d);
        }
#endif
        ws->stp = min(ws->stp, stpmax);
        opl_csrch_start(&ws->lnsrch, *f, ws->gd, ws->stp,
                        SFTOL(ws), SGTOL(ws), SXTOL(ws),
                        0, stpmax);
        if (TASK(ws) != OPL_TASK_FG) {
            /* Some error occurred (task and context should have been set so as
               to describe the problem). */
            break;
        }
        ws->searching = OPL_TRUE;
        /* Compute the new iterate. */
        return next_step(ws, x);

    case OPL_TASK_ERROR:
        /* Nothing to do then. */
        break;

    default:
        /* Probably an error. */
        return failure(ws, OPL_CORRUPTED, "corrupted workspace");
    }

    return TASK(ws);

# undef Y
# undef S
# undef SUCCESS
# undef FAILURE

}

/*
 * Compute new search direction H(k).d(k) based on the two-loop recursion
 * L-BFGS formula (operation is done in-place).  H(k) is the limited memory
 * BFGS approximation of the inverse Hessian, d(k) has been initialized with is
 * the (projected) gradient at k-th step.  H(k) is approximated by using the M
 * last pairs (s, y) where:
 *
 *   s(j) = x(j+1) - x(j)
 *   y(j) = g(j+1) - g(j)
 *
 * or:
 *
 *   s(j) = x(j) - x(j+1)
 *   y(j) = g(j) - g(j+1)
 *
 * Indeed, the way the differences are computed does not matter (as far as all
 * differences are computed similarly), because: it has no influence on RHO(j),
 * and on dot products between S(j) and Y(j); it changes the sign of ALPHA(j)
 * and BETA but not that of ALPHA(j) - BETA times S(j) or Y(j).
 *
 * The two-loop recursion algorithm writes:
 *
 *   1- start with current gradient:
 *        d := +/- g(k)
 *
 *   2- for j = k-1, ..., k-m
 *        rho(j) = s(j)'.y(j)           (save rho(j))
 *        alpha(j) = (s(j)'.d)/rho(j)   (save alpha(j))
 *        d := d - alpha(j) y(j)
 *
 *   3- apply approximation of inverse k-th Hessian:
 *        d := H0(k).d
 *      for instance:
 *        d := (rho(k-1)/(y(k-1)'.y(k-1))) d
 *
 *   4- for j = k-m, ..., k-1
 *        d := d + (alpha(j) - (y(j)'.d)/rho(j)) s(j)
 */
static void compute_direction(
    opl_vmlmb_workspace_t* ws,
    double d[],
    const int isfree[],
    const double h[])
{
    double gamma = (isfree == NULL ? ws->gamma : 0);
    long m = ws->m;
    long mp = ws->mp;
    long n = ws->n;
    long off = ws->mark + m + 1;
    double* rho = ws->rho;
    double* alpha = ws->alpha;
    double** S_ = ws->S;
    double** Y_ = ws->Y;

# define S(j) S_[j]
# define Y(j) Y_[j]

    /* First loop of the L-BFGS recursion. */
    for (long k = 1; k <= mp; ++k) {
        long j = (off - k)%m;
        if (isfree != NULL) {
            rho[j] = opl_ddot_free(n, S(j), Y(j), isfree);
        }
        if (rho[j] > 0) {
            alpha[j] = opl_ddot(n, S(j), d)/rho[j];
            opl_daxpy_free(n, -alpha[j], Y(j), d, isfree);
            if (gamma <= 0) {
                double yty = opl_ddot_free(n, Y(j), Y(j), isfree);
                if (yty > 0) {
                    gamma = rho[j]/yty;
                }
            }
        }
    }

    /* Apply initial approximation of the inverse Hessian. */
    if (h != NULL) {
        /* Apply diagonal preconditioner. */
        for (long i = 0; i < n; ++i) {
            d[i] *= h[i];
        }
    } else if (gamma > 0) {
        /* Apply initial H (i.e. just scale D) and perform the second stage of
           the 2-loop recursion. */
        opl_dscal(n, gamma, d);
    } else {
        /* All correction pairs are invalid: manage to restart the L-BFGS
           recursion.  Note that D is left unchanged in that case. */
        ++ws->restarts;
        ws->mp = 0;
        return;
    }

    /* Second loop of the L-BFGS recursion. */
    for (long k = mp; k >= 1; --k) {
        long j = (off - k)%m;
        if (rho[j] > 0) {
            double beta = opl_ddot(n, Y(j), d)/rho[j];
            opl_daxpy_free(n, alpha[j] - beta, S(j), d, isfree);
        }
    }

# undef S
# undef Y

}

static opl_task_t next_step(
    opl_vmlmb_workspace_t* ws,
    double x[])
{
    double stp = ws->stp;
    const double* d = ws->d;
    const double* x0 = ws->S[ws->mark];
    long n = ws->n;

    /* Compute the new iterate. */
    for (long i = 0; i < n; ++i) {
        x[i] = x0[i] - stp*d[i];
    }

    /* Require the caller to compute the function and its gradient. */
    return success(ws, OPL_TASK_FG, "compute f(x) and g(x)");
}

static int check_free(
    opl_vmlmb_workspace_t* ws,
    int isfree[],
    const double h[])
{
    if (h != NULL) {
        const double zero = 0;
        long n = ws->n;
        if (isfree != NULL) {
            /* fix ISFREE array */
            for (long i = 0; i < n; ++i) {
                if (isfree[i] && h[i] <= zero) {
                    isfree[i] = 0;
                }
            }
        } else {
            /* check that H is positive definite */
            for (long i = 0; i < n; ++i) {
                if (h[i] <= zero) {
                    failure(
                        ws, OPL_NOT_POSITIVE_DEFINITE,
                        "initial inverse Hessian is not positive definite");
                    return -1;
                }
            }
        }
    }
    return 0;
}

/*---------------------------------------------------------------------------*/
/* MEMORY MANAGEMENT FOR VMLMB */

/*
 * A monolithic workspace is organized in memory as follows (starting from the
 * base address):
 *
 *   base:           workspace structure
 *                   padding for proper alignment
 *   base + offset1: 2 arrays of M pointers (for S and Y)
 *                   padding for proper alignment
 *   base + offset2: 2 arrays of M reals (for ALPHA and RHO)
 *                   1 array of N reals (for D, the search direction)
 *                   2*M arrays of N reals (for S and Y)
 */

static const size_t offset1 = OPL_ROUND_UP(sizeof(opl_vmlmb_workspace_t),
                                           sizeof(void*));

static size_t workspace_offset2(
    size_t m)
{
    return OPL_ROUND_UP(offset1 + 2*m*sizeof(void*), sizeof(double));
}

static size_t number_of_reals(
    size_t n,
    size_t m)
{
    return (2*m*(n + 1) + n);
}

size_t opl_vmlmb_monolithic_workspace_size(
    long n,
    long m)
{
    if (n <= 0 || m <= 0) {
        errno = EINVAL;
        return 0;
    }
    return workspace_offset2(m) + number_of_reals(n, m)*sizeof(double);
}

opl_vmlmb_workspace_t* opl_vmlmb_monolithic_workspace_init(
    void* buf,
    long n,
    long m)
{
    /* Check arguments. */
    if (buf == NULL) {
        if (errno != ENOMEM) {
            errno = EFAULT;
        }
        return NULL;
    }
    if (n <= 0 || m <= 0) {
        errno = EINVAL;
        return NULL;
    }

    /* Compute offsets and size.  Clear buffer. */
    size_t offset2 = workspace_offset2(m);
    size_t size = offset2 + number_of_reals(n, m)*sizeof(double);
    memset(buf, 0, size);

    /* Instanciate workspace. */
    opl_vmlmb_workspace_t* ws = (opl_vmlmb_workspace_t*)buf;
    ws->m = m;
    ws->n = n;
    ws->S = (double**)((byte_t*)ws + offset1);
    ws->Y = ws->S + m;
    ws->alpha = (double*)((byte_t*)ws + offset2);
    ws->rho = ws->alpha + m;
    ws->d = ws->rho + m;
    double* arr = ws->d;
    for (long k = 0; k < m; ++k) {
        ws->S[k] = (arr += n);
        ws->Y[k] = (arr += n);
    }
    opl_vmlmb_restart(ws);
    return opl_vmlmb_set_defaults(ws);
}

opl_vmlmb_workspace_t* opl_vmlmb_set_defaults(
    opl_vmlmb_workspace_t* ws)
{
    SFTOL(ws) = DEFAULT_SFTOL;
    SGTOL(ws) = DEFAULT_SGTOL;
    SXTOL(ws) = DEFAULT_SXTOL;
    ws->fatol = DEFAULT_FATOL;
    ws->frtol = DEFAULT_FRTOL;
    ws->delta = DEFAULT_DELTA;
    ws->epsilon = DEFAULT_EPSILON;
    opl_vmlmb_set_fmin(ws, NaN);
    return ws;
}

static void free_split_workspace(
    void* ptr)
{
    opl_vmlmb_workspace_t* ws = (opl_vmlmb_workspace_t*)ptr;
    long m = ws->m;
    double* tmp;

    if ((tmp = ws->d) != NULL) {
        ws->d = NULL;
        free(tmp);
    }
    for (long k = 0; k < m; ++k) {
        if ((tmp = ws->S[k]) != NULL) {
            ws->S[k] = NULL;
            free(tmp);
        }
        if ((tmp = ws->Y[k]) != NULL) {
            ws->Y[k] = NULL;
            free(tmp);
        }
    }
    free(ws);
}

opl_vmlmb_workspace_t* opl_vmlmb_create(
    long n,
    long m)
{
    if (n <= 0 || m <= 0) {
        errno = EINVAL;
        return 0;
    }
    if (m*n <= 10000) {
        /* For small problems, the workspace is allocated as a single block of
           memory. */
        size_t size = opl_vmlmb_monolithic_workspace_size(n, m);
        opl_vmlmb_workspace_t* ws = opl_vmlmb_monolithic_workspace_init(
            malloc(size), n, m);
        if (ws == NULL) {
            return NULL;
        }
        ws->free = free;
        return ws;
    } else {
        /* For larger problems, the philosophy is to store the workspace
           structure and related small arrays in a single block of memory,
           while larger arrays (d, S and Y) are stored separately. */
        size_t offset2 = workspace_offset2(m);
        size_t size = offset2 + 2*m*sizeof(double);
        opl_vmlmb_workspace_t* ws = (opl_vmlmb_workspace_t*)malloc(size);
        if (ws == NULL) {
            return NULL;
        }
        memset(ws, 0, size);
        ws->free = free_split_workspace;
        ws->m = m;
        ws->n = n;
        ws->S = (double**)((byte_t*)ws + offset1);
        ws->Y = ws->S + m;
        ws->alpha = (double*)((byte_t*)ws + offset2);
        ws->rho = ws->alpha + m;
        if ((ws->d = NEW(double, n)) == NULL) {
        destroy:
            opl_vmlmb_destroy(ws);
            return NULL;
        }
        for (long k = 0; k < m; ++k) {
            if ((ws->S[k] = NEW(double, n)) == NULL ||
                (ws->Y[k] = NEW(double, n)) == NULL) {
                goto destroy;
            }
        }
        opl_vmlmb_restart(ws);
        return opl_vmlmb_set_defaults(ws);
    }
}

void opl_vmlmb_destroy(
    opl_vmlmb_workspace_t* ws)
{
    if (ws != NULL) {
        if (ws->free == NULL) {
            fprintf(stderr, "*** ERROR *** %s\n",
                    "attempt to destroy a foreign monolithic workspace, "
                    "read the documentation!");
        } else {
            ws->free(ws);
        }
    }
}
