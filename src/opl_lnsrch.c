/*
 * opl_lnsrch.c --
 *
 * Line search routines for OptimPackLegacy library.
 *
 *-----------------------------------------------------------------------------
 *
 * This file is part of OptimPackLegacy
 * <https://github.com/emmt/OptimPackLegacy>.
 *
 * Copyright (c) 2003-2019, Éric Thiébaut.
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
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include "opl_private.h"

/* Constants. */
static const double zero = 0.0;
static const double xtrapl = 1.1;
static const double xtrapu = 4.0;

/*---------------------------------------------------------------------------*/

const size_t OPL_CSRCH_WORKSPACE_SIZE =
  OPL_ROUND_UP(sizeof(opl_csrch_workspace_t), sizeof(double));

size_t
opl_csrch_get_workspace_size()
{
  return OPL_CSRCH_WORKSPACE_SIZE;
}

opl_csrch_workspace_t*
opl_csrch_create_workspace()
{
  opl_csrch_workspace_t* ws = OPL_NEW(opl_csrch_workspace_t, 1);
  if (ws != NULL) {
    memset(ws, 0, sizeof(opl_csrch_workspace_t));
    opl_initialize_context(&ws->context);
  }
  return ws;
}

void
opl_csrch_destroy_workspace(opl_csrch_workspace_t* ws)
{
  if (ws != NULL) {
    free((void*)ws);
  }
}

opl_task_t
opl_csrch_get_task(opl_csrch_workspace_t* ws)
{
  return (ws != NULL ? ws->task : OPL_TASK_ERROR);
}

opl_status_t
opl_csrch_get_status(opl_csrch_workspace_t* ws)
{
  return (ws != NULL ? ws->context.status : OPL_ILLEGAL_ADDRESS);
}

const char*
opl_csrch_get_reason(opl_csrch_workspace_t* ws)
{
  return (ws != NULL ? ws->context.message : "illegal workspace address");
}

static opl_status_t
report(opl_csrch_workspace_t* ws,
       opl_task_t task, opl_status_t status, const char* mesg)
{
  ws->task = task;
  return opl_set_context((opl_context_t*)ws, status, mesg, OPL_PERMANENT);
}

opl_status_t
opl_csrch_start(opl_csrch_workspace_t* ws, double f, double g, double stp,
                double ftol, double gtol, double xtol,
                double stpmin, double stpmax)
{
  /* Minimal checking. */
  if (ws == NULL) {
    errno = EFAULT;
    return OPL_ILLEGAL_ADDRESS;
  }

  /* Check the input arguments for errors.
     Exit if there are errors on input. */
#define ERROR(I,S) report(ws, OPL_TASK_ERROR, I, "opl_csrch_start: " S)
    if (stpmax < stpmin) {
      return ERROR(OPL_STPMAX_LT_STPMIN, "STPMAX < STPMIN");
    }
    if (stpmin < zero) {
      return ERROR(OPL_STPMIN_LT_ZERO, "STPMIN < 0");
    }
    if (xtol < zero) {
      return ERROR(OPL_XTOL_LT_ZERO, "XTOL < 0");
    }
    if (ftol <= zero) {
      return ERROR(OPL_FTOL_LE_ZERO, "FTOL <= 0");
    }
    if (gtol <= zero) {
      return ERROR(OPL_GTOL_LE_ZERO, "GTOL <= 0");
    }
    if (g >= zero) {
      return ERROR(OPL_NOT_A_DESCENT, "initial G >= 0");
    }
    if (stp > stpmax) {
      return ERROR(OPL_STP_GT_STPMAX, "STP > STPMAX");
    }
    if (stp < stpmin) {
      return ERROR(OPL_STP_LT_STPMIN, "STP < STPMIN");
    }
#undef ERROR

    /* Initialize local variables.
       The variables STX, FX, GX contain the values of the step,
       function, and derivative at the best step.
       The variables STY, FY, GY contain the value of the step,
       function, and derivative at STY.
       The variables STP, F, G contain the values of the step,
       function, and derivative at STP. */
    ws->ftol    = ftol;
    ws->gtol    = gtol;
    ws->xtol    = xtol;
    ws->stpmin  = stpmin;
    ws->stpmax  = stpmax;
    ws->finit   = f;
    ws->ginit   = g;
    ws->stx     = zero;
    ws->fx      = f;
    ws->gx      = g;
    ws->sty     = zero;
    ws->fy      = f;
    ws->gy      = g;
    ws->stmin   = zero;
    ws->stmax   = stp + stp*xtrapu;
    ws->width   = stpmax - stpmin;
    ws->width1  = 2.0*(stpmax - stpmin);
    ws->brackt  = 0;
    ws->stage   = 1;
    ws->task    = OPL_TASK_FG;
    return opl_success((opl_context_t*)ws);
}

opl_status_t
opl_csrch_iterate(opl_csrch_workspace_t* ws,
                  double f, double g, double *stp_ptr)
{
  /* Define macros to make the code simpler and easier to read. */
# define ftol    (ws->ftol)
# define gtol    (ws->gtol)
# define xtol    (ws->xtol)
# define stpmin  (ws->stpmin)
# define stpmax  (ws->stpmax)
# define finit   (ws->finit)
# define ginit   (ws->ginit)
# define stx     (ws->stx)
# define fx      (ws->fx)
# define gx      (ws->gx)
# define sty     (ws->sty)
# define fy      (ws->fy)
# define gy      (ws->gy)
# define stmin   (ws->stmin)
# define stmax   (ws->stmax)
# define width   (ws->width)
# define width1  (ws->width1)
# define task    (ws->task)
# define stage   (ws->stage)
# define brackt  (ws->brackt)

  /* Local variables. */
  double ftest, gtest, stp;
  opl_status_t status;

  /* Minimal checking. */
  if (ws == NULL) {
    errno = EFAULT;
    return OPL_ILLEGAL_ADDRESS;
  }

  /* Initialize local variables. */
  stp = *stp_ptr;

  /* If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the algorithm enters
     the second stage. */
  gtest = ftol*ginit;
  ftest = finit + stp*gtest;
  if (stage == 1 && f <= ftest && g >= zero) {
    stage = 2;
  }

  /* Test for termination: convergence or warnings. */
  if (f <= ftest && fabs(g) <= gtol*(-ginit)) {
    return report(ws, OPL_TASK_CONV,
                  OPL_SUCCESS, "strong Wolfe conditions both satisfied");
  }
  if (stp == stpmin && (f > ftest || g >= gtest)) {
    return report(ws, OPL_TASK_WARN,
                  OPL_STP_EQ_STPMIN, "step at lower bound");
  }
  if (stp == stpmax && f <= ftest && g <= gtest) {
    return report(ws, OPL_TASK_WARN,
                  OPL_STP_EQ_STPMAX, "step at upper bound");
  }
  if (brackt && stmax - stmin <= xtol*stmax) {
    return report(ws, OPL_TASK_WARN,
                  OPL_XTOL_TEST_SATISFIED, "XTOL test satisfied");
  }
  if (brackt && (stp <= stmin || stp >= stmax)) {
    return report(ws, OPL_TASK_WARN,
                  OPL_ROUNDING_ERROR, "rounding errors prevent progress");
  }

  /* A modified function is used to predict the step during the first
     stage if a lower function value has been obtained but the decrease
     is not sufficient. */

  if (stage == 1 && f <= fx && f > ftest) {

    /* Define the modified function and derivative values. */
    double fm  = f  - stp*gtest;
    double fxm = fx - stx*gtest;
    double fym = fy - sty*gtest;
    double gm  = g  - gtest;
    double gxm = gx - gtest;
    double gym = gy - gtest;

    /* Call opl_cstep to update STX, STY, and to compute the new step. */
    status = opl_cstep((opl_context_t*)ws, &brackt, stmin, stmax,
                       &stx, &fxm, &gxm,
                       &sty, &fym, &gym,
                       &stp,  fm,   gm);
    if (status != OPL_SUCCESS) {
      task = OPL_TASK_ERROR;
      return status;
    }

    /* Reset the function and derivative values for F. */
    fx = fxm + stx*gtest;
    fy = fym + sty*gtest;
    gx = gxm + gtest;
    gy = gym + gtest;

  } else {

    /* Call opl_cstep to update STX, STY, and to compute the new step. */
    status = opl_cstep((opl_context_t*)ws, &brackt, stmin, stmax,
                       &stx, &fx, &gx,
                       &sty, &fy, &gy,
                       &stp,  f,   g);
    if (status != OPL_SUCCESS) {
      task = OPL_TASK_ERROR;
      return status;
    }

  }

  /* Decide if a bisection step is needed. */
  if (brackt) {
    double wcur = fabs(sty - stx);
    if (wcur >= 0.66*width1) stp = stx + 0.5*(sty - stx);
    width1 = width;
    width = wcur;
  }

  /* Set the minimum and maximum steps allowed for STP. */
  if (brackt) {
    if (stx <= sty) {
      stmin = stx;
      stmax = sty;
    } else {
      stmin = sty;
      stmax = stx;
    }
  } else {
    stmin = stp + xtrapl*(stp - stx);
    stmax = stp + xtrapu*(stp - stx);
  }

  /* Force the step to be within the bounds STPMAX and STPMIN. */
  if (stp > stpmax) stp = stpmax;
  if (stp < stpmin) stp = stpmin;

  /* If further progress is not possible, let STP be the best
     point obtained during the search. */
  if (brackt && (stp <= stmin || stp >= stmax ||
                 stmax - stmin <= xtol*stmax)) {
    stp = stx;
  }

  /* Save local variables. */
  *stp_ptr = stp;

  /* Obtain another function and derivative. */
  return report(ws, OPL_TASK_FG,
                OPL_SUCCESS, "compute f(x) and g(x)");

  /* Undefine macros. */
# undef ftol
# undef gtol
# undef xtol
# undef stpmin
# undef stpmax
# undef finit
# undef ginit
# undef stx
# undef fx
# undef gx
# undef sty
# undef fy
# undef gy
# undef stmin
# undef stmax
# undef width
# undef width1
# undef task
# undef stage
# undef brackt
}

/*---------------------------------------------------------------------------*/
/* SAFEGUARDED CUBIC STEP */

opl_status_t
opl_cstep(opl_context_t* ctx, opl_boolean_t *brackt,
          double stpmin, double stpmax,
          double *stx_ptr, double *fx_ptr, double *dx_ptr,
          double *sty_ptr, double *fy_ptr, double *dy_ptr,
          double *stp_ptr, double  fp,     double  dp)
{
  /* Get values of input/output variables. */
  double stx = *stx_ptr, fx = *fx_ptr, dx = *dx_ptr;
  double sty = *sty_ptr, fy = *fy_ptr, dy = *dy_ptr;
  double stp = *stp_ptr;

  /* Local variables. */
  double gamma, theta, p, q, r, s, t;
  double stpc; /* cubic step */
  double stpq; /* quadratic step */
  double sgnd, stpf;

  /* Check the input parameters for errors. */
# define REPORT(status, mesg)  opl_set_context(ctx, status,             \
                                               "opl_cstep: " mesg,      \
                                               OPL_PERMANENT)
  if (*brackt && (stx < sty ? (stp <= stx || stp >= sty)
                  : (stp >= stx || stp <= sty))) {
    return REPORT(OPL_OUT_OF_BOUNDS, "STP outside bracket (STX,STY)");
  } else if (dx*(stp - stx) >= zero) {
    return REPORT(OPL_NOT_A_DESCENT, "descent condition violated");
  } else if (stpmax < stpmin) {
    return REPORT(OPL_STPMAX_LT_STPMIN, "STPMAX < STPMIN");
    return 0;
  }
# undef REPORT

  /* Determine if the derivatives have opposite sign. */
  sgnd = (dx/fabs(dx))*dp;

  if (fp > fx) {
    /* First case.  A higher function value.  The minimum is bracketed.  If the
       cubic step is closer to STX than the quadratic step, the cubic step is
       taken, otherwise the average of the cubic and quadratic steps is
       taken. */
    theta = 3.0*(fx - fp)/(stp - stx) + dx + dp;
    s = fabs(theta);
    if (s < (t = fabs(dx))) s = t;
    if (s < (t = fabs(dp))) s = t;
    t = theta/s;
    gamma = s*sqrt(t*t - (dx/s)*(dp/s));
    if (stp < stx) gamma = -gamma;
    p = (gamma - dx) + theta;
    q = ((gamma - dx) + gamma) + dp;
    r = p/q;
    stpc = stx + r*(stp - stx);
    stpq = stx + ((dx/((fx - fp)/(stp - stx) + dx))/2.0)*(stp - stx);
    if (fabs(stpc - stx) < fabs(stpq - stx)) {
      stpf = stpc;
    } else {
      stpf = stpc + (stpq - stpc)/2.0;
    }
    *brackt = 1;
  } else if (sgnd < zero) {
    /* Second case.  A lower function value and derivatives of opposite sign.
       The minimum is bracketed.  If the cubic step is farther from STP than
       the secant (quadratic) step, the cubic step is taken, otherwise the
       secant step is taken. */
    theta = 3.0*(fx - fp)/(stp - stx) + dx + dp;
    s = fabs(theta);
    if (s < (t = fabs(dx))) s = t;
    if (s < (t = fabs(dp))) s = t;
    t = theta/s;
    gamma = s*sqrt(t*t - (dx/s)*(dp/s));
    if (stp > stx) gamma = -gamma;
    p = (gamma - dp) + theta;
    q = ((gamma - dp) + gamma) + dx;
    r = p/q;
    stpc = stp + r*(stx - stp);
    stpq = stp + (dp/(dp - dx))*(stx - stp);
    if (fabs(stpc - stp) > fabs(stpq - stp)) {
      stpf = stpc;
    } else {
      stpf = stpq;
    }
    *brackt = 1;
  } else if (fabs(dp) < fabs(dx)) {
    /* Third case.  A lower function value, derivatives of the same sign, and
       the magnitude of the derivative decreases.  The cubic step is computed
       only if the cubic tends to infinity in the direction of the step or if
       the minimum of the cubic is beyond STP.  Otherwise the cubic step is
       defined to be the secant step. */
    theta = 3.0*(fx - fp)/(stp - stx) + dx + dp;
    s = fabs(theta);
    if (s < (t = fabs(dx))) s = t;
    if (s < (t = fabs(dp))) s = t;
    /* The case GAMMA = 0 only arises if the cubic does not tend to infinity in
       the direction of the step. */
    t = theta/s;
    t = t*t - (dx/s)*(dp/s);
    gamma = (t > zero ? s*sqrt(t) : zero);
    if (stp > stx) gamma = -gamma;
    p = (gamma - dp) + theta;
    /*q = ((gamma - dp) + gamma) + dx;*/
    q = (gamma + (dx - dp)) + gamma;
    r = p/q;
    if (r < zero && gamma != zero) {
      stpc = stp + r*(stx - stp);
    } else if (stp > stx) {
      stpc = stpmax;
    } else {
      stpc = stpmin;
    }
    stpq = stp + (dp/(dp - dx))*(stx - stp);
    if (*brackt) {
      /* A minimizer has been bracketed.  If the cubic step is closer to STP
	 than the secant step, the cubic step is taken, otherwise the secant
	 step is taken. */
      stpf = fabs(stpc - stp) < fabs(stpq - stp) ? stpc : stpq;
      t = stp + 0.66*(sty - stp);
      if (stp > stx ? stpf > t : stpf < t) stpf = t;
    } else {
      /* A minimizer has not been bracketed. If the cubic step is farther from
	 stp than the secant step, the cubic step is taken, otherwise the
	 secant step is taken. */
      stpf = fabs(stpc - stp) > fabs(stpq - stp) ? stpc : stpq;
      if (stpf > stpmax) stpf = stpmax;
      if (stpf < stpmin) stpf = stpmin;
    }
  } else {
    /* Fourth case.  A lower function value, derivatives of the same sign, and
       the magnitude of the derivative does not decrease.  If the minimum is
       not bracketed, the step is either STPMIN or STPMAX, otherwise the cubic
       step is taken. */
    if (*brackt) {
      theta = 3.0*(fp - fy)/(sty - stp) + dy + dp;
      s = fabs(theta);
      if (s < (t = fabs(dy))) s = t;
      if (s < (t = fabs(dp))) s = t;
      t = theta/s;
      gamma = s*sqrt(t*t - (dy/s)*(dp/s));
      if (stp > sty) gamma = -gamma;
      p = (gamma - dp) + theta;
      q = ((gamma - dp) + gamma) + dy;
      r = p/q;
      stpc = stp + r*(sty - stp);
      stpf = stpc;
    } else if (stp > stx) {
      stpf = stpmax;
    } else {
      stpf = stpmin;
    }
  }

  /* Update the interval which contains a minimizer and guess for next step. */
  if (fp > fx) {
    *sty_ptr = stp;
    *fy_ptr = fp;
    *dy_ptr = dp;
  } else {
    if (sgnd < zero) {
      *sty_ptr = stx;
      *fy_ptr = fx;
      *dy_ptr = dx;
    }
    *stx_ptr = stp;
    *fx_ptr = fp;
    *dx_ptr = dp;
  }
  *stp_ptr = stpf;
  return opl_success(ctx);
}

/*---------------------------------------------------------------------------*/
