/*
 * op_vmlmb --
 *
 *	Variable Metric Limited Memory with Bound constraints for
 *	OptimPack library.
 *
 *-----------------------------------------------------------------------------
 *
 *      Copyright (c) 2003, Eric THIEBAUT.
 *
 *	This file is part of OptimPack.
 *
 *	OptimPack is  free software; you can redistribute  it and/or modify
 *	it under the  terms of the GNU General  Public License as published
 *	by the Free  Software Foundation; either version 2  of the License,
 *	or (at your option) any later version.
 *
 *	OptimPack is  distributed in the hope  that it will  be useful, but
 *	WITHOUT  ANY  WARRANTY;  without   even  the  implied  warranty  of
 *	MERCHANTABILITY or  FITNESS FOR A PARTICULAR PURPOSE.   See the GNU
 *	General Public License for more details.
 *
 *	You should have  received a copy of the  GNU General Public License
 *	along with OptimPack (file  "LICENSE" in the top source directory);
 *	if  not, write  to the  Free Software  Foundation, Inc.,  59 Temple
 *	Place, Suite 330, Boston, MA 02111-1307 USA
 *
 *-----------------------------------------------------------------------------
 *
 *	$Id$
 *	$Log$
 *
 *-----------------------------------------------------------------------------
 */

#include <math.h>
#include <string.h>
#include <stdio.h>
#include "optimpack.h"

#define op_dcopy(n, x, y) memcpy((y), (x), (n)*sizeof(double))

/*---------------------------------------------------------------------------*/
/* LIMITED MEMORY VARIABLE METRIC METHOD (BFGS)
   WITH/WITHOUT BOUND CONSTRAINTS */

/* Indices 0-1 in ISAVE are reserved for op_cshrc. */
#define INDEX_OF_TASK   2
#define INDEX_OF_STAGE  3
#define INDEX_OF_M      4
#define INDEX_OF_N      5
#define INDEX_OF_ITER   6
#define INDEX_OF_MARK   7
#define INDEX_OF_MP     8

/* Indices 0-11 in ISAVE are reserved for op_cshrc. */
#define INDEX_OF_SFTOL  12
#define INDEX_OF_SGTOL  13
#define INDEX_OF_SXTOL  14
#define INDEX_OF_FRTOL  15
#define INDEX_OF_FATOL  16
#define INDEX_OF_FMIN   17
#define INDEX_OF_F0     18
#define INDEX_OF_GD     19
#define INDEX_OF_GD0    20
#define INDEX_OF_STP    21
#define INDEX_OF_STPMIN 22
#define INDEX_OF_STPMAX 23
#define INDEX_OF_WORK   24 /* must be the last one */
/* DSAVE must have at least: 24 + 2×M + N + 2×N×M elements */

#define SET_TASK(val, str) (op_mcopy("op_vmlmb_next: " str, csave), task=(val))

int op_vmlmb_first(op_integer_t n, op_integer_t m,
		   double fmin, double fatol, double frtol,
		   double sftol, double sgtol, double sxtol,
		   char csave[], op_integer_t isave[], double dsave[])
{
  int task = OP_TASK_FG;
  if (n <= 0) return SET_TASK(OP_TASK_ERROR, "N <= 0");
  if (m <= 0) return SET_TASK(OP_TASK_ERROR, "M <= 0");
  if (fatol < 0.0) return SET_TASK(OP_TASK_ERROR, "FATOL < 0");
  if (frtol < 0.0) return SET_TASK(OP_TASK_ERROR, "FRTOL < 0");
  if (sxtol <= 0.0) return SET_TASK(OP_TASK_ERROR, "SXTOL <= 0");
  if (sxtol >= 1.0) return SET_TASK(OP_TASK_ERROR, "SXTOL >= 1");
  if (sftol <= 0.0) return SET_TASK(OP_TASK_ERROR, "SFTOL <= 0");
  if (sftol >= 1.0) return SET_TASK(OP_TASK_ERROR, "SFTOL >= 1");
  if (sgtol <= 0.0) return SET_TASK(OP_TASK_ERROR, "SGTOL <= 0");
  if (sgtol >= 1.0) return SET_TASK(OP_TASK_ERROR, "SGTOL >= 1");
  if (sftol >= sgtol) return SET_TASK(OP_TASK_ERROR, "SFTOL >= SGTOL");

  isave[INDEX_OF_TASK]  = OP_TASK_FG;
  isave[INDEX_OF_STAGE] = 0;
  isave[INDEX_OF_M]     = m;
  isave[INDEX_OF_N]     = n;
  dsave[INDEX_OF_FMIN]  = fmin;
  dsave[INDEX_OF_FRTOL] = frtol;
  dsave[INDEX_OF_FATOL] = fatol;
  dsave[INDEX_OF_SFTOL] = sftol;
  dsave[INDEX_OF_SGTOL] = sgtol;
  dsave[INDEX_OF_SXTOL] = sxtol;
  return isave[INDEX_OF_TASK];
}
#undef SET_TASK

#define S(K)     &s[(K)*n]
#define Y(K)     &y[(K)*n]
#define SET_TASK(val, str) (op_mcopy("op_vmlmb_next: " str, csave), task=(val))

int op_vmlmb_next(double x[], double *f, double g[], const int active[],
		  char csave[], op_integer_t isave[], double dsave[])
{
  const double zero = 0.0;

  /* Get local variables.  (Note: depending on the value of STAGE, it is
     possible to restore only _some_ of theses variables, but this would
     result in a less readable code with negligible speedup gain.) */
  int          task  = isave[INDEX_OF_TASK];
  int          stage = isave[INDEX_OF_STAGE];
  op_integer_t m     = isave[INDEX_OF_M];
  op_integer_t n     = isave[INDEX_OF_N];
  op_integer_t iter  = isave[INDEX_OF_ITER];
  op_integer_t mark  = isave[INDEX_OF_MARK];
  op_integer_t mp    = isave[INDEX_OF_MP];

  double sftol  = dsave[INDEX_OF_SFTOL];
  double sgtol  = dsave[INDEX_OF_SGTOL];
  double sxtol  = dsave[INDEX_OF_SXTOL];
  double fmin   = dsave[INDEX_OF_FMIN];
  double frtol  = dsave[INDEX_OF_FRTOL];
  double fatol  = dsave[INDEX_OF_FATOL];
  double f0     = dsave[INDEX_OF_F0];
  double gd     = dsave[INDEX_OF_GD];
  double gd0    = dsave[INDEX_OF_GD0];
  double stp    = dsave[INDEX_OF_STP];
  double stpmin = dsave[INDEX_OF_STPMIN];
  double stpmax = dsave[INDEX_OF_STPMAX];

  double *alpha = dsave + INDEX_OF_WORK;
  double *rho = alpha + m;
  double *x0  = rho + m;
  double *s   = x0 + n;
  double *y   = s + n*m;
  double *ptr;
  double q;
  int info;
  op_integer_t i, j, k;

  /*
   * Differences with original FORTRAN version:
   *
   *   (5) The effective step is used instead of the search direction
   *       times the step size (useful, e.g., when parameter constraints
   *       are imposed by projections).
   *
   *   (6) It is possible to use an active subset of parameters, e.g. to apply
   *       bound constraints.
   *
   *   (7) Discard pairs for which RHO = S'.Y <= 0 (i.e. H must be positive
   *       definite).
   */


  /* STAGE is: 0 - first entry
   *           1 - start line search
   *           2 - line search in progress
   *           3 - line search converged */


  /* S(MARK) = -D(k) where D(k) is the k-th search direction */

  if (stage == 0) {
    /* First search direction is the (normalized) steepest descent. */
    if (*f <= fmin) {
      SET_TASK(OP_TASK_ERROR, "initial F <= FMIN");
      goto done;
    }
    iter = 0;  /* number of successful iterations */
  restart:
    mark = 0;  /* index of current direction */
    mp = 0;    /* number of saved directions */
  steepest:
    op_dcopy_active(n, g, S(mark), active); /* steepest ascent */
    q = op_dnrm2(n, S(mark)); /* same as |g| but account for active */
    if (q > zero) {
      op_dscal(n, 1.0/q, S(mark));
      gd = -q;
    } else {
      SET_TASK(OP_TASK_CONV, "local minimum found");
      goto done;
    }
    stage = 1; /* set STAGE to initialize the line search */

  } else {

    if (stage == 3) {
      /* Previous step was successful.  Compute new search direction
       * H(k).g(k) based on the two-loop recursion L-BFGS formula.  H(k) is
       * the limited memory BFGS approximation of the inverse Hessian, g(k)
       * is the gradient at k-th step.  H(k) is approximated by using the M
       * last pairs (s, y) where:
       *
       *   s(j) = x(j+1) - x(j)
       *   y(j) = g(j+1) - g(j)
       *
       * The two-loop recursion algorithm reads:
       *
       *   1- start with current gradient:
       *        v := g(k)
       *
       *   2- for j = k-1, ..., k-m
       *        rho(j) = s(j)'.y(j)           (save rho(j))
       *        alpha(j) = (s(j)'.v)/rho(j)   (save alpha(j))
       *        v := v - alpha(j) y(j)
       *
       *   3- apply approximation of inverse k-th Hesian:
       *        v := H0(k).v
       *      for instance:
       *        v := (rho(k-1)/(y(k-1)'.y(k-1))) v
       *
       *   4- for j = k-m, ..., k-1
       *        v := v + (alpha(j) - (y(j)'.v)/rho(j)) s(j)
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
      double *v = x0;
      double gamma = zero;
      op_integer_t mm = mark + m;
      op_dcopy_active(n, g, v, active);
      for (k=0 ; k<mp ; ++k) {
	j = (mm - k)%m;
	if (active) rho[j] = op_ddot_active(n, S(j), Y(j), active);
	if (rho[j] > zero) {
	  alpha[j] = op_ddot(n, S(j), v)/rho[j];
	  op_daxpy_active(n, -alpha[j], Y(j), v, active);
	  if (! gamma) gamma = rho[j]/op_ddot_active(n, Y(j), Y(j), active);
	}
      }
      if (gamma) {
	/* Apply initial H (i.e. just scale V) and perform the second stage of
	   the 2-loop recursion. */
	op_dscal(n, gamma, v);
      } else {
	/* All correction pairs are invalid: fallback to use the steepest
	   descent direction. */
	fprintf(stderr, "WARNING: %s\n",
		"no valid correction pair (use steepest descent direction)");
	goto steepest;
      }
      for (k=mp-1 ; k>=0 ; --k) {
	j = (mm - k)%m;
	if (rho[j] <= zero) continue;
	op_daxpy_active(n, alpha[j] -  op_ddot(n, Y(j), v)/rho[j],
			S(j), v, active);
      }

      /* Save search direction. */
      mark = (mark + 1)%m;
      op_dcopy(n, v, S(mark));

      /* Set STAGE to initialize the line search. */
      stage = 1;
    }

    /* Compute derivative with respect to step size STP. */
    gd = -op_ddot(n, g, S(mark));
  }

  if (stage == 1) {
    /* Set variables so as to initialize the line search subroutine. */
    if (gd >= zero) {
      /* L-BFGS recursion yields a search direction which is not a descent.
         We therefore restart the algorithm with steepest descent. */
      fprintf(stderr, "WARNING: %s\n",
	      "not a descent direction (algorithm restarted)");
      goto restart;
    }
    f0 = *f;
    gd0 = gd;
    stpmin = zero;
    stpmax = (fmin - f0)/(sgtol*gd0);
    stp = OP_MIN(1.0, stpmax);
    op_dcopy(n, x, x0);      /* save parameters */
    op_dcopy(n, g, Y(mark)); /* save gradient */
    stage = 2;
    task = OP_TASK_START; /* set TASK so as to start line search */
  } else {
    task = OP_TASK_FG; /* set TASK so as to continue line search */
  }

  if (stage == 2) {
    /* Determine the line search parameter. */
    if (*f < fmin) {
      SET_TASK(OP_TASK_WARN, "F < FMIN");
    } else {
      /* Call line search iterator. */
      info = op_csrch(*f, gd, &stp, sftol, sgtol, sxtol, stpmin, stpmax, &task,
		      csave, isave, dsave);
      if (info == 1) {
        /* Compute the new iterate. */
	for (ptr=S(mark), i=0 ; i<n ; ++i) x[i] = x0[i] - stp*ptr[i];
      } else if (info == 2 || info == 5) {
	/* Line search has converged. */
	++iter;
        if (mp < m) ++mp;
	stage = 3;

        /* Compute the step and gradient change.  Note: we always compute
	   the effective step (parameter difference) to account for bound
	   constraints and, at least, numerical rounding errors. */
	for (ptr=Y(mark), i=0 ; i<n ; ++i) ptr[i] -= g[i];
	for (ptr=S(mark), i=0 ; i<n ; ++i) ptr[i] = x0[i] - x[i];
	if (! active) rho[mark] = op_ddot(n, Y(mark), S(mark));

	/* Test for global convergence otherwise set TASK to signal a new
	   iterate.  Set STAGE to compute a new search direction. */
        if (op_noneof(n, S(mark))) {
	  SET_TASK(OP_TASK_WARN, "no parameter change");
        } else if (op_noneof(n, Y(mark))) {
	  SET_TASK(OP_TASK_WARN, "no gradient change");
        } else {
	  double change1 = fabs(*f - f0);
	  double change2 = fabs(stp*gd0);
          double change = OP_MAX(change1, change2);
          if (change <= frtol*fabs(f0)) {
            SET_TASK(OP_TASK_CONV, "FRTOL test satisfied");
          } else if (change <= fatol) {
            SET_TASK(OP_TASK_CONV, "FATOL test satisfied");
          } else {
            SET_TASK(OP_TASK_NEWX,
		     "new improved solution available for inspection");
          }
        }
      } else if (info >= 3) {
        /* INFO >= 3 means that line search could not converge.
	   Restore solution at start of line search. */
        op_dcopy(n, x0, x);
        op_dcopy(n, Y(mark), g);
        *f = f0;
      }
    }
  }

  /* Save local variables (but constant ones) and return TASK. */
 done:
#warning "task not needed"
  isave[INDEX_OF_TASK]   = task;
  isave[INDEX_OF_STAGE]  = stage;
  isave[INDEX_OF_ITER]   = iter;
  isave[INDEX_OF_MARK]   = mark;
  isave[INDEX_OF_MP]     = mp;
#if 0
  dsave[INDEX_OF_SFTOL]  = sftol; /* constant */
  dsave[INDEX_OF_SGTOL]  = sgtol; /* constant */
  dsave[INDEX_OF_SXTOL]  = sxtol; /* constant */
  dsave[INDEX_OF_FRTOL]  = frtol; /* constant */
  dsave[INDEX_OF_FATOL]  = fatol; /* constant */
  dsave[INDEX_OF_FMIN]   = fmin;  /* constant */
#endif
  dsave[INDEX_OF_F0]     = f0;
  dsave[INDEX_OF_GD]     = gd;
  dsave[INDEX_OF_GD0]    = gd0;
  dsave[INDEX_OF_STP]    = stp;
  dsave[INDEX_OF_STPMIN] = stpmin;
  dsave[INDEX_OF_STPMAX] = stpmax;
  return task;
}

#undef SET_TASK
#undef Y
#undef S

double op_vmlmb_set_fmin(const char csave[], const op_integer_t isave[],
			 double dsave[], double new_value)
{
  double old_value = dsave[INDEX_OF_FMIN];
  dsave[INDEX_OF_FMIN] = new_value;
  return old_value;
}

#define GET_FUNC(type, foo, FOO, save)				\
type OP_CONCAT(op_vmlmb_get_,foo)(const char csave[],		\
				  const op_integer_t isave[],	\
				  const double dsave[])		\
{ return dsave[OP_CONCAT(INDEX_OF_,FOO)]; }
GET_FUNC(double, sftol, SFTOL, dsave)
GET_FUNC(double, sgtol, SGTOL, dsave)
GET_FUNC(double, sxtol, SXTOL, dsave)
GET_FUNC(double, frtol, FRTOL, dsave)
GET_FUNC(double, fatol, FATOL, dsave)
GET_FUNC(double, fmin,  FMIN, dsave)
GET_FUNC(double, step,  STP, dsave)
GET_FUNC(op_integer_t, iter,  ITER, isave)
#undef GET_FUNC

/*---------------------------------------------------------------------------*/
