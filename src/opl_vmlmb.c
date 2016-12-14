/*
 * opl_vmlmb --
 *
 * Variable Metric Limited Memory with Bound constraints for OptimPackLegacy
 * library.
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (c) 2003-2009, 2016 Éric Thiébaut.
 *
 * This file is part of OptimPack <https://github.com/emmt/OptimPackLegacy>.
 *
 * OptimPack is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option)
 * any later version.
 *
 * OptimPack is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * OptimPack (file "LICENSE" in the top source directory); if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307 USA
 *
 *-----------------------------------------------------------------------------
 */

#include <math.h>
#include <string.h>
#include <stdio.h>
#include "optimpacklegacy.h"

#define opl_dcopy(n, x, y) memcpy((y), (x), (n)*sizeof(double))

#define STPMAX 1e20

const double zero = 0.0;
const double one = 1.0;

static double min(double a, double b) {
  return (a <= b ? a : b);
}

static double max(double a, double b) {
  return (a >= b ? a : b);
}

/*---------------------------------------------------------------------------*/
/* LIMITED MEMORY VARIABLE METRIC METHOD (BFGS) WITH/WITHOUT BOUND
   CONSTRAINTS */

/* Indices 0-1 in ISAVE are reserved for opl_cshrc. */
#define INDEX_OF_TASK       2
#define INDEX_OF_M          3
#define INDEX_OF_N          4
#define INDEX_OF_ITER       5
#define INDEX_OF_MARK       6
#define INDEX_OF_MP         7
#define INDEX_OF_FLAGS      8
#define INDEX_OF_NEVALS     9
#define INDEX_OF_NRESTARTS 10
#if (OPL_VMLMB_ISAVE_NUMBER != INDEX_OF_NRESTARTS + 1)
# error bad ISAVE settings
#endif

/* Indices 0-11 in ISAVE are reserved for opl_cshrc. */
#define INDEX_OF_SFTOL    12
#define INDEX_OF_SGTOL    13
#define INDEX_OF_SXTOL    14
#define INDEX_OF_FRTOL    15
#define INDEX_OF_FATOL    16
#define INDEX_OF_FMIN     17
#define INDEX_OF_F0       18
#define INDEX_OF_GD       19
#define INDEX_OF_G0D      20
#define INDEX_OF_STP      21
#define INDEX_OF_STPMIN   22
#define INDEX_OF_STPMAX   23
#define INDEX_OF_DELTA    24
#define INDEX_OF_EPSILON  25
#define INDEX_OF_GNORM    26
#define INDEX_OF_G0NORM   27
#define INDEX_OF_WORK     28 /* must be the last one */
#if (OPL_VMLMB_DSAVE_NUMBER(0,0) != INDEX_OF_WORK)
# error bad DSAVE settings
#endif

#define FLAG(n)                          (1 << (n))
#define FLAG_IS_SET(value, flag)         (((value) & (flag)) != 0)
#define FLAG_FMIN                        FLAG(0)

#define ERROR(str)   SET_TASK(OPL_TASK_ERROR, str)

#define VMLMB_SETUP "opl_vmlmb_setup: "
#define VMLMB_NEXT  "opl_vmlmb_next: "

int
opl_vmlmb_setup(opl_integer_t n, opl_integer_t m,
                double fatol, double frtol,
                double sftol, double sgtol, double sxtol,
                double delta, double epsilon,
                char csave[], opl_integer_t isave[], double dsave[])
{
#define SET_TASK(val, str) (opl_mcopy(VMLMB_SETUP str, csave), task = (val))
  int task = OPL_TASK_FG;
  if (n <= 0) return ERROR("N <= 0");
  if (m <= 0) return ERROR("M <= 0");
  if (fatol < 0.0) return ERROR("FATOL < 0");
  if (frtol < 0.0) return ERROR("FRTOL < 0");
  if (sxtol <= 0.0) return ERROR("SXTOL <= 0");
  if (sxtol >= 1.0) return ERROR("SXTOL >= 1");
  if (sftol <= 0.0) return ERROR("SFTOL <= 0");
  if (sftol >= 1.0) return ERROR("SFTOL >= 1");
  if (sgtol <= 0.0) return ERROR("SGTOL <= 0");
  if (sgtol >= 1.0) return ERROR("SGTOL >= 1");
  if (sftol >= sgtol) return ERROR("SFTOL >= SGTOL");
  if (delta < 0.0) return ERROR("DELTA < 0");
  if (epsilon < 0.0) return ERROR("EPSILON < 0");

  csave[0]                  = '\0';

  isave[INDEX_OF_TASK]      =  OPL_TASK_FG;
  isave[INDEX_OF_M]         =  m;
  isave[INDEX_OF_N]         =  n;
  isave[INDEX_OF_ITER]      =  0;
  isave[INDEX_OF_MARK]      = -1;
  isave[INDEX_OF_MP]        =  0;
  isave[INDEX_OF_FLAGS]     =  0;
  isave[INDEX_OF_NEVALS]    =  0;
  isave[INDEX_OF_NRESTARTS] =  0;

  dsave[INDEX_OF_SFTOL]     =  sftol;
  dsave[INDEX_OF_SGTOL]     =  sgtol;
  dsave[INDEX_OF_SXTOL]     =  sxtol;
  dsave[INDEX_OF_FRTOL]     =  frtol;
  dsave[INDEX_OF_FATOL]     =  fatol;
  dsave[INDEX_OF_FMIN]      =  0.0;
  dsave[INDEX_OF_F0]        =  0.0;
  dsave[INDEX_OF_GD]        =  0.0;
  dsave[INDEX_OF_G0D]       =  0.0;
  dsave[INDEX_OF_STP]       =  0.0;
  dsave[INDEX_OF_STPMIN]    =  0.0;
  dsave[INDEX_OF_STPMAX]    =  0.0;
  dsave[INDEX_OF_DELTA]     =  delta;
  dsave[INDEX_OF_EPSILON]   =  epsilon;
  dsave[INDEX_OF_GNORM]     = -1.0;
  dsave[INDEX_OF_G0NORM]    = -1.0;

  return isave[INDEX_OF_TASK];
#undef SET_TASK
}

int
opl_vmlmb_restart(char csave[], opl_integer_t isave[], double dsave[])
{
  csave[0]                  = '\0';

  isave[INDEX_OF_TASK]      =  OPL_TASK_FG;
  isave[INDEX_OF_ITER]      =  0;
  isave[INDEX_OF_MARK]      = -1;
  isave[INDEX_OF_MP]        =  0;
  isave[INDEX_OF_NEVALS]    =  0;
  isave[INDEX_OF_NRESTARTS] =  0;

  dsave[INDEX_OF_F0]        =  0.0;
  dsave[INDEX_OF_GD]        =  0.0;
  dsave[INDEX_OF_G0D]       =  0.0;
  dsave[INDEX_OF_STP]       =  0.0;
  dsave[INDEX_OF_STPMIN]    =  0.0;
  dsave[INDEX_OF_STPMAX]    =  0.0;
  dsave[INDEX_OF_GNORM]     = -1.0;
  dsave[INDEX_OF_G0NORM]    = -1.0;

  return isave[INDEX_OF_TASK];
}

int
opl_vmlmb_restore(double x[], double *f, double g[],
                  char csave[], opl_integer_t isave[], double dsave[])
{
  int task = isave[INDEX_OF_TASK];
  if (task == OPL_TASK_FG) {
    opl_integer_t m = isave[INDEX_OF_M];
    opl_integer_t n = isave[INDEX_OF_N];
    opl_integer_t mark = isave[INDEX_OF_MARK];
    const double *s = dsave + INDEX_OF_WORK + 2*m + n;
    const double *y = s + n*m;
    *f = dsave[INDEX_OF_F0];
    dsave[INDEX_OF_GNORM] = dsave[INDEX_OF_G0NORM];
    opl_dcopy(n, &s[mark*n], x);
    opl_dcopy(n, &y[mark*n], g);
    task = OPL_TASK_NEWX;
    isave[INDEX_OF_TASK] = task;
  }
  return task;
}

/*---------------------------------------------------------------------------*/

/* Check H and ISFREE arrays, possibly fix ISFREE, returns non-zero
   in case of error. */
static int
check_free(opl_integer_t n, opl_logical_t isfree[], const double h[],
           int *task, char csave[]);

/* Compute next step and set task. */
static void
next_step(opl_integer_t n, double x[], const double x0[], double stp,
          const double d[], int *task, char csave[]);

/* A bunch of macros to simplify the code. */
#define SET_TASK(val, str) (opl_mcopy(VMLMB_NEXT str, csave), task=(val))
#define S(K)               &s[(K)*n]
#define Y(K)               &y[(K)*n]

int
opl_vmlmb_next(double x[], double *f, double g[],
               opl_logical_t isfree[], const double h[],
               char csave[], opl_integer_t isave[], double dsave[])
{
  /* Get local variables. */
  int           task      = isave[INDEX_OF_TASK];
  opl_integer_t m         = isave[INDEX_OF_M];
  opl_integer_t n         = isave[INDEX_OF_N];
  opl_integer_t iter      = isave[INDEX_OF_ITER];
  opl_integer_t mark      = isave[INDEX_OF_MARK];
  opl_integer_t mp        = isave[INDEX_OF_MP];
  opl_integer_t flags     = isave[INDEX_OF_FLAGS];
  opl_integer_t nevals    = isave[INDEX_OF_NEVALS];
  opl_integer_t nrestarts = isave[INDEX_OF_NRESTARTS];
  int           have_fmin = ((flags & FLAG_FMIN) != 0);

  double sftol    = dsave[INDEX_OF_SFTOL];
  double sgtol    = dsave[INDEX_OF_SGTOL];
  double sxtol    = dsave[INDEX_OF_SXTOL];
  double fmin     = dsave[INDEX_OF_FMIN];
  double frtol    = dsave[INDEX_OF_FRTOL];
  double fatol    = dsave[INDEX_OF_FATOL];
  double f0       = dsave[INDEX_OF_F0];
  double gd       = dsave[INDEX_OF_GD];
  double g0d      = dsave[INDEX_OF_G0D];
  double stp      = dsave[INDEX_OF_STP];
  double stpmin   = dsave[INDEX_OF_STPMIN];
  double stpmax   = dsave[INDEX_OF_STPMAX];
  double delta    = dsave[INDEX_OF_DELTA];
  double epsilon  = dsave[INDEX_OF_EPSILON];
  double gnorm    = dsave[INDEX_OF_GNORM];
  double g0norm   = dsave[INDEX_OF_G0NORM];

  double *alpha = dsave + INDEX_OF_WORK;
  double *rho = alpha + m;
  double *d   = rho + m; /* anti-search direction */
  double *s   = d + n;
  double *y   = s + n*m;
  double *Sm, *Ym;
  int info, descent;
  opl_integer_t i, j, k;

  switch (task) {

  case OPL_TASK_FG:
    /* Caller has perfomed a new evaluation of the function and its gradient. */
    ++nevals;
    if (have_fmin && *f <= fmin) {
      ERROR("initial F <= FMIN");
      break;
    }
    if (nevals > 1) {
      /* Compute directional derivative for the line search. */
      gd = -opl_ddot(n, g, d);

      /* Call line search iterator to check for line search convergence. */
      info = opl_csrch(*f, gd, &stp, sftol, sgtol, sxtol, stpmin, stpmax,
                       &task, csave, isave, dsave);
      if (info == 2 || info == 5) {
        /* Line search has converged. */
        ++iter;
      } else {
        /* Line search has not converged. */
        if (info == 1) {
          /* Compute the new iterate. */
          next_step(n, x, S(mark), stp, d, &task, csave);
        } else if (info != 2 && info != 5) {
          /* Error or warning in line search (TASK and CSAVE should be set so
             as to describe the problem).  Restore solution at start of line
             search. */
          opl_dcopy(n, S(mark), x);
          opl_dcopy(n, Y(mark), g);
          gnorm = g0norm;
          *f = f0;
        }
        break;
      }
    }

    /* Request the caller to compute the set of unbinded variables. */
    SET_TASK(OPL_TASK_FREEVARS, "determine free variables");
    break;

  case OPL_TASK_FREEVARS:
    /* Copy the (projected) gradient into D and compute the norm of the
       (projected) gradient. */
    if (check_free(n, isfree, h, &task, csave) != 0) {
      break;
    }
    opl_dcopy_free(n, g, d, isfree);
    gnorm = opl_dnrm2(n, d);
    if (gnorm == zero) {
      SET_TASK(OPL_TASK_CONV, "local minimum found");
      break;
    }
    if (nevals > 1) {
      /* Compute the step and gradient change (i.e., update L-BFGS data).
         Note: we always compute the effective step (parameter difference) to
         account for bound constraints and, at least, numerical rounding or
         truncation errors. */
      for (Ym = Y(mark), i = 0; i < n; ++i) {
        Ym[i] -= g[i];
      }
      for (Sm = S(mark), i = 0; i < n; ++i) {
        Sm[i] -= x[i];
      }
      if (isfree == NULL) {
        rho[mark] = opl_ddot(n, Y(mark), S(mark));
      }
      if (mp < m) {
        ++mp;
      }

      /* Test for global convergence. */
      if (opl_noneof(n, S(mark))) {
        SET_TASK(OPL_TASK_WARN, "no parameter change");
      } else if (opl_noneof(n, Y(mark))) {
        SET_TASK(OPL_TASK_WARN, "no gradient change");
      } else {
        double change = max(fabs(*f - f0), fabs(stp*g0d));
        if (change <= frtol*fabs(f0)) {
          SET_TASK(OPL_TASK_CONV, "FRTOL test satisfied");
        } else if (change <= fatol) {
          SET_TASK(OPL_TASK_CONV, "FATOL test satisfied");
        }
      }
    }
    if (task == OPL_TASK_FREEVARS) {
      /* Set task to signal a new iterate. */
      SET_TASK(OPL_TASK_NEWX, "new improved solution available for inspection");
    }
    break;

  case OPL_TASK_NEWX:
  case OPL_TASK_CONV:
  case OPL_TASK_WARN:
    /* Compute a new search direction.  D already contains the (projected) gradient. */
    if (mp > 0) {
      /* Compute new search direction H(k).g(k)
       * based on the two-loop recursion L-BFGS formula.  H(k) is the limited
       * memory BFGS approximation of the inverse Hessian, g(k) is the
       * gradient at k-th step.  H(k) is approximated by using the M last
       * pairs (s, y) where:
       *
       *   s(j) = x(j+1) - x(j)
       *   y(j) = g(j+1) - g(j)
       *
       * The two-loop recursion algorithm writes:
       *
       *   1- start with current gradient:
       *        d := g(k)
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
       *
       * Note: in fact, this program saves -S and -Y, but by looking at
       * above relations it is clear that this change of sign:
       *  - has no influence on RHO, and dot products between S and Y;
       *  - changes the sign of ALPHA but not that of ALPHA(j) times
       *    S(j) or Y(j);
       * the two-loop recursion therefore involves the same equations
       * with:     s(j) = x(j+1) - x(j)  and  y(j) = g(j+1) - g(j)
       * or with:  s(j) = x(j) - x(j+1)  and  y(j) = g(j) - g(j+1).
       */
      double gamma = zero;
      opl_integer_t off = mark + m + 1;

      for (k = 1; k <= mp; ++k) {
        j = (off - k)%m;
        if (isfree != NULL) {
          rho[j] = opl_ddot_free(n, S(j), Y(j), isfree);
        }
        if (rho[j] > zero) {
          alpha[j] = opl_ddot(n, S(j), d)/rho[j];
          opl_daxpy_free(n, -alpha[j], Y(j), d, isfree);
          if (gamma <= zero) {
            double yty = opl_ddot_free(n, Y(j), Y(j), isfree);
            if (yty > zero) {
              gamma = rho[j]/yty;
            }
          }
        }
      }
      if (h != NULL) {
        /* Apply diagonal preconditioner. */
        for (i = 0; i < n; ++i) {
          d[i] *= h[i];
        }
      } else if (gamma > zero) {
        /* Apply initial H (i.e. just scale D) and perform the second stage of
           the 2-loop recursion. */
        opl_dscal(n, gamma, d);
      } else {
        /* All correction pairs are invalid: manage to restart the L-BFGS
           recursion.  Note that D has not changed so it is the (projected)
           gradient. */
        ++nrestarts;
        mp = 0;
      }
      for (k = mp; k >= 1; --k) {
        j = (off - k)%m;
        if (rho[j] > zero) {
          opl_daxpy_free(n, alpha[j] - opl_ddot(n, Y(j), d)/rho[j],
                         S(j), d, isfree);
        }
      }

      if (mp > 0) {
        /* Set initial step size and compute dot product of gradient and search
           direction.  If L-BFGS recursion filed to produce a sufficient descent
           search direction, we restart the algorithm with steepest descent. */
        stp = 1.0;
        gd = -opl_ddot(n, g, d);
        if (epsilon > zero) {
          descent = (gd <= -epsilon*gnorm*opl_dnrm2(n, d));
        } else {
          descent = (gd < zero);
        }
        if (! descent) {
          /* Manage to restart the L-BFGS with steepest descent. */
          ++nrestarts;
          mp = 0;
          opl_dcopy_free(n, g, d, isfree);
        }
      }
    }

    if (mp == 0) {
      /* Compute initial search direction (or after a restart).  D is already
         the (projected) gradient. */
      if (h != NULL) {
        /* Use diagonal preconditioner to compute initial search direction. */
        for (i = 0; i < n; ++i) {
          d[i] *= h[i];
        }
        stp = 1.0;
        gd = -opl_ddot(n, g, d);
        if (gd >= zero) {
          ERROR("preconditioner is not positive definite");
          break;
        }
      } else {
        /* No preconditioning, use a small step along the steepest descent. */
        if (delta > zero) {
          stp = (opl_dnrm2(n, x)/gnorm)*delta;
        } else {
          stp = zero;
        }
        if (stp <= zero) {
          /* The following step length is quite arbitrary but to avoid this
             would require to know the typical size of the variables. */
          stp = one/gnorm;
        }
        gd = -gnorm*gnorm;
      }
    }

    /* Advance the mark and save point at start of line search.  Note that this
       point has forcibly been projected so it is feasible. */
    mark = (mark + 1)%m;
    f0 = *f;
    g0d = gd;
    g0norm = gnorm;
    opl_dcopy(n, x, S(mark)); /* save parameters X0 */
    opl_dcopy(n, g, Y(mark)); /* save gradient G0 */

    /* Set step bounds and task to start line search. */
    stpmin = zero;
#if 1
    stpmax = STPMAX;
#else
    if (have_fmin) {
      stpmax = (fmin - f0)/(sgtol*g0d);
    } else {
      double temp = fabs(f0);
      if (temp < 1.0) {
	temp = 1.0;
      }
      stpmax = temp/(sgtol*g0d);
    }
#endif
    stp = min(stp, stpmax);
    task = OPL_TASK_START;
    info = opl_csrch(*f, gd, &stp, sftol, sgtol, sxtol, stpmin, stpmax,
                     &task, csave, isave, dsave);
    if (info == 1) {
      /* Compute the new iterate (otherwise, it must be an error and there is
         nothing to do). */
      next_step(n, x, S(mark), stp, d, &task, csave);
    }
    break;

  case OPL_TASK_ERROR:
    /* Nothing to do then. */
    break;

  default:
    /* Probably an error. */
    SET_TASK(OPL_TASK_ERROR, "corrupted workspace");
    break;
  }

  /* Save local variables (but constant ones) and return TASK. */
  isave[INDEX_OF_TASK]      = task;
  isave[INDEX_OF_ITER]      = iter;
  isave[INDEX_OF_MARK]      = mark;
  isave[INDEX_OF_MP]        = mp;
#if 0
  isave[INDEX_OF_FLAGS]     = flags; /* constant */
#endif
  isave[INDEX_OF_NEVALS]    = nevals;
  isave[INDEX_OF_NRESTARTS] = nrestarts;

#if 0
  dsave[INDEX_OF_SFTOL]     = sftol; /* constant */
  dsave[INDEX_OF_SGTOL]     = sgtol; /* constant */
  dsave[INDEX_OF_SXTOL]     = sxtol; /* constant */
  dsave[INDEX_OF_FRTOL]     = frtol; /* constant */
  dsave[INDEX_OF_FATOL]     = fatol; /* constant */
  dsave[INDEX_OF_FMIN]      = fmin;  /* constant */
#endif
  dsave[INDEX_OF_F0]        = f0;
  dsave[INDEX_OF_GD]        = gd;
  dsave[INDEX_OF_G0D]       = g0d;
  dsave[INDEX_OF_STP]       = stp;
  dsave[INDEX_OF_STPMIN]    = stpmin;
  dsave[INDEX_OF_STPMAX]    = stpmax;
#if 0
  dsave[INDEX_OF_DELTA]     = delta;   /* constant */
  dsave[INDEX_OF_EPSILON]   = epsilon; /* constant */
#endif
  dsave[INDEX_OF_GNORM]    = gnorm;
  dsave[INDEX_OF_G0NORM]   = g0norm;
  return task;
}

#undef SET_TASK
#undef Y
#undef S

/*---------------------------------------------------------------------------*/

static void
next_step(opl_integer_t n, double x[], const double x0[], double stp,
          const double d[], int* task, char csave[])
{
  opl_integer_t i;

  /* Compute the new iterate. */
  for (i = 0; i < n; ++i) {
    x[i] = x0[i] - stp*d[i];
  }

  *task = OPL_TASK_FG;
  opl_mcopy(VMLMB_NEXT "compute f(x) and g(x)", csave);
}

static int
check_free(opl_integer_t n, opl_logical_t isfree[],
           const double h[], int* task, char csave[])
{
  if (h != NULL) {
    const double zero = 0.0;
    opl_integer_t i;
    if (isfree != NULL) {
      /* fix ISFREE array */
      for (i = 0; i < n; ++i) {
	if (isfree[i] && h[i] <= zero) {
	  isfree[i] = 0;
	}
      }
    } else {
      /* check that H is positive definite */
      for (i = 0; i < n; ++i) {
	if (h[i] <= zero) {
	  opl_mcopy("opl_vmlmb_next: H is not positive definite", csave);
	  *task = OPL_TASK_ERROR;
	  return -1;
	}
      }
    }
  }
  return 0;
}

/*---------------------------------------------------------------------------*/

#define EMIT_CODE(name, NAME)                                   \
  int opl_vmlmb_set_##name(const char csave[],                  \
                           opl_integer_t isave[],               \
                           double dsave[],                      \
                           double new_value,                    \
                           double *old_value)                   \
  {                                                             \
    int test = ((isave[INDEX_OF_FLAGS] & FLAG_##NAME) != 0);    \
    if (test && old_value != (double *)NULL) {                  \
      *old_value = dsave[INDEX_OF_##NAME];                      \
    }                                                           \
    dsave[INDEX_OF_##NAME] = new_value;                         \
    isave[INDEX_OF_FLAGS] |= FLAG_##NAME;                       \
    return test;                                                \
  }                                                             \
  int opl_vmlmb_get_##name(const char csave[],                  \
                           const opl_integer_t isave[],         \
                           const double dsave[],                \
                           double *ptr)                         \
  {                                                             \
    int test = ((isave[INDEX_OF_FLAGS] & FLAG_##NAME) != 0);    \
    if (test && ptr != (double *)NULL) {                        \
      *ptr = dsave[INDEX_OF_##NAME];                            \
    }                                                           \
    return test;                                                \
  }
EMIT_CODE(fmin, FMIN)
#undef EMIT_CODE

#define EMIT_CODE(type, foo, FOO, save)                                 \
  type OPL_CONCAT(opl_vmlmb_get_,foo)(const char csave[],               \
                                      const opl_integer_t isave[],      \
                                      const double dsave[])             \
  {                                                                     \
    return dsave[OPL_CONCAT(INDEX_OF_,FOO)];                            \
  }
EMIT_CODE(double, sftol,      SFTOL,    dsave)
EMIT_CODE(double, sgtol,      SGTOL,    dsave)
EMIT_CODE(double, sxtol,      SXTOL,    dsave)
EMIT_CODE(double, frtol,      FRTOL,    dsave)
EMIT_CODE(double, fatol,      FATOL,    dsave)
EMIT_CODE(double, step,       STP,      dsave)
EMIT_CODE(double, delta,      DELTA,    dsave)
EMIT_CODE(double, epsilon,    EPSILON,  dsave)
EMIT_CODE(double, gnorm,     GNORM,   dsave)
EMIT_CODE(opl_integer_t, iter,      ITER,      isave)
EMIT_CODE(opl_integer_t, nevals,    NEVALS,    isave)
EMIT_CODE(opl_integer_t, nrestarts, NRESTARTS, isave)
#undef EMIT_CODE

/*---------------------------------------------------------------------------*/
