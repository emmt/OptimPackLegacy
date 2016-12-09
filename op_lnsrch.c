/*
 * op_lnsrch.h --
 *
 *	Line search routines for OptimPack library.
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
 *	$Id: op_lnsrch.c,v 1.1 2003/03/12 15:11:52 eric Exp $
 *	$Log: op_lnsrch.c,v $
 *	Revision 1.1  2003/03/12 15:11:52  eric
 *	Initial revision
 *
 *
 *-----------------------------------------------------------------------------
 */

#include <math.h>
#include <string.h>
#include "optimpack.h"

/*---------------------------------------------------------------------------*/

/* Indices of variables saved in ISAVE array: */
#define INDEX_OF_BRACKT    0
#define INDEX_OF_STAGE     1

/* Indices of variables saved in DSAVE array: */
#define INDEX_OF_FINIT     0
#define INDEX_OF_GINIT     1
#define INDEX_OF_STX       2
#define INDEX_OF_FX        3
#define INDEX_OF_GX        4
#define INDEX_OF_STY       5
#define INDEX_OF_FY        6
#define INDEX_OF_GY        7
#define INDEX_OF_STMIN     8
#define INDEX_OF_STMAX     9
#define INDEX_OF_WIDTH    10
#define INDEX_OF_WIDTH1   11

#define SET_TASK(val, str) *task = (val); op_mcopy("op_csrch: " str, csave)

int op_csrch(double f, double g, double *stp_ptr,
	     double ftol, double gtol, double xtol,
	     double stpmin, double stpmax, int *task,
	     char csave[], op_integer_t isave[], double dsave[])
{
  const double zero = 0.0;
  const double xtrapl = 1.1;
  const double xtrapu = 4.0;
  int brackt, stage, info;
  double finit, ginit, ftest, gtest;
  double fx, gx, stx;
  double fy, gy, sty;
  double width, width1;
  double stmin, stmax;
  double stp = *stp_ptr;

  if (*task == OP_TASK_START) {
    /* Check the input arguments for errors.
       Exit if there are errors on input. */
#define RETURN_ERROR(I,S) SET_TASK(OP_TASK_ERROR, S); return (I)
    if (stpmax < stpmin) { RETURN_ERROR( 0, "STPMAX < STPMIN"); }
    if (stpmin < zero)   { RETURN_ERROR(-3,"STPMIN < 0"); }
    if (xtol < zero)     { RETURN_ERROR(-4, "XTOL < 0"); }
    if (ftol <= zero)    { RETURN_ERROR(-5, "FTOL <= 0"); }
    if (gtol <= zero)    { RETURN_ERROR(-6, "GTOL <= 0"); }
    if (g >= zero)       { RETURN_ERROR(-7, "initial G >= 0"); }
    if (stp > stpmax)    { RETURN_ERROR(-8, "STP > STPMAX"); }
    if (stp < stpmin)    { RETURN_ERROR(-9, "STP < STPMIN"); }
#undef RETURN_ERROR

    /* Initialize local variables.
       The variables STX, FX, GX contain the values of the step,
       function, and derivative at the best step.
       The variables STY, FY, GY contain the value of the step,
       function, and derivative at STY.
       The variables STP, F, G contain the values of the step,
       function, and derivative at STP. */
    isave[INDEX_OF_BRACKT] = 0;
    isave[INDEX_OF_STAGE]  = 1;
    dsave[INDEX_OF_FINIT]  = f;
    dsave[INDEX_OF_GINIT]  = g;
    dsave[INDEX_OF_STX]    = zero;
    dsave[INDEX_OF_FX]     = f;
    dsave[INDEX_OF_GX]     = g;
    dsave[INDEX_OF_STY]    = zero;
    dsave[INDEX_OF_FY]     = f;
    dsave[INDEX_OF_GY]     = g;
    dsave[INDEX_OF_STMIN]  = zero;
    dsave[INDEX_OF_STMAX]  = stp + stp*xtrapu;
    dsave[INDEX_OF_WIDTH]  = stpmax - stpmin;
    dsave[INDEX_OF_WIDTH1] = 2.0*(stpmax - stpmin);
    *task = OP_TASK_FG;
    return 1;
  }

  /* Restore local variables. */
  brackt = isave[INDEX_OF_BRACKT];
  stage  = isave[INDEX_OF_STAGE];
  finit  = dsave[INDEX_OF_FINIT];
  ginit  = dsave[INDEX_OF_GINIT];
  stx    = dsave[INDEX_OF_STX];
  fx     = dsave[INDEX_OF_FX];
  gx     = dsave[INDEX_OF_GX];
  sty    = dsave[INDEX_OF_STY];
  fy     = dsave[INDEX_OF_FY];
  gy     = dsave[INDEX_OF_GY];
  stmin  = dsave[INDEX_OF_STMIN];
  stmax  = dsave[INDEX_OF_STMAX];
  width  = dsave[INDEX_OF_WIDTH];
  width1 = dsave[INDEX_OF_WIDTH1];

  /* If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
     algorithm enters the second stage. */
  gtest = ftol*ginit;
  ftest = finit + stp*gtest;
  if (stage == 1 && f <= ftest && g >= zero) stage = 2;

  /* Test for termination: convergence or warnings. */
  if (f <= ftest && fabs(g) <= gtol*(-ginit)) {
    /* Strong Wolfe conditions both satisfied. */
    SET_TASK(OP_TASK_CONV, "convergence of line search");
    info = 2;
  } else if (stp == stpmin && (f > ftest || g >= gtest)) {
    SET_TASK(OP_TASK_WARN, "STP = STPMIN");
    info = 3;
  } else if (stp == stpmax && f <= ftest && g <= gtest) {
    SET_TASK(OP_TASK_WARN, "STP = STPMAX");
    info = 4;
  } else if (brackt && stmax - stmin <= xtol*stmax) {
    SET_TASK(OP_TASK_WARN, "XTOL test satisfied");
    info = 5;
  } else if (brackt && (stp <= stmin || stp >= stmax)) {
    SET_TASK(OP_TASK_WARN, "rounding errors prevent progress");
    info = 6;
  } else {

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

      /* Call op_cstep to update STX, STY, and to compute the new step. */
      info = op_cstep(&stx, &fxm, &gxm,
		      &sty, &fym, &gym,
		      &stp,  fm,   gm,
		      &brackt, stmin, stmax, csave);
      if (info <= 0) {
	*task = OP_TASK_ERROR;
	return info;
      }

      /* Reset the function and derivative values for F. */
      fx = fxm + stx*gtest;
      fy = fym + sty*gtest;
      gx = gxm + gtest;
      gy = gym + gtest;

    } else {

      /* Call op_cstep to update STX, STY, and to compute the new step. */
      info = op_cstep(&stx, &fx, &gx,
		      &sty, &fy, &gy,
		      &stp,  f,   g,
		      &brackt, stmin, stmax, csave);
      if (info <= 0) {
	*task = OP_TASK_ERROR;
	return info;
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
		   stmax - stmin <= xtol*stmax)) stp = stx;

    /* Obtain another function and derivative. */
    if (csave) csave[0] = 0;
    *task = OP_TASK_FG;
    info = 1;
  }

  /* Save local variables (only the ones that may have changed). */
  *stp_ptr = stp;
  isave[INDEX_OF_BRACKT] = brackt;
  isave[INDEX_OF_STAGE]  = stage;
  dsave[INDEX_OF_STX]    = stx;
  dsave[INDEX_OF_FX]     = fx;
  dsave[INDEX_OF_GX]     = gx;
  dsave[INDEX_OF_STY]    = sty;
  dsave[INDEX_OF_FY]     = fy;
  dsave[INDEX_OF_GY]     = gy;
  dsave[INDEX_OF_STMIN]  = stmin;
  dsave[INDEX_OF_STMAX]  = stmax;
  dsave[INDEX_OF_WIDTH]  = width;
  dsave[INDEX_OF_WIDTH1] = width1;
  return info;
}

#undef stp
#undef SET_TASK

/*---------------------------------------------------------------------------*/

int op_cstep(double *stx_ptr, double *fx_ptr, double *dx_ptr,
	     double *sty_ptr, double *fy_ptr, double *dy_ptr,
	     double *stp_ptr, double  fp,     double  dp,
	     int *brackt, double stpmin, double stpmax,
	     char csave[])
{
  /* Constants. */
  const double zero = 0.0;

  /* Get values of input/output variables. */
  double stx = *stx_ptr, fx = *fx_ptr, dx = *dx_ptr;
  double sty = *sty_ptr, fy = *fy_ptr, dy = *dy_ptr;
  double stp = *stp_ptr;

  /* Local variables. */
  double gamma, theta, p, q, r, s, t;
  double stpc; /* cubic step */
  double stpq; /* quadratic step */
  double sgnd, stpf;
  int info;

  /* Check the input parameters for errors. */
  if (*brackt && (stx < sty ? (stp <= stx || stp >= sty)
		            : (stp >= stx || stp <= sty))) {
    op_mcopy("op_cstep: STP outside bracket (STX,STY)", csave);
    return -2;
  } else if (dx*(stp - stx) >= zero) {
    op_mcopy("op_cstep: descent condition violated", csave);
    return -1;
  } else if (stpmax < stpmin) {
    op_mcopy("op_cstep: STPMAX < STPMIN", csave);
    return 0;
  }

  /* Determine if the derivatives have opposite sign. */
  sgnd = (dx/fabs(dx))*dp;

  if (fp > fx) {
    /* First case.  A higher function value.  The minimum is bracketed.  If
       the cubic step is closer to STX than the quadratic step, the cubic
       step is taken, otherwise the average of the cubic and quadratic
       steps is taken. */
    info = 1;
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
    /* Second case.  A lower function value and derivatives of opposite
       sign.  The minimum is bracketed.  If the cubic step is farther from
       STP than the secant (quadratic) step, the cubic step is taken,
       otherwise the secant step is taken. */
    info = 2;
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
    /* Third case.  A lower function value, derivatives of the same sign,
       and the magnitude of the derivative decreases.  The cubic step is
       computed only if the cubic tends to infinity in the direction of the
       step or if the minimum of the cubic is beyond STP.  Otherwise the
       cubic step is defined to be the secant step. */
    info = 3;
    theta = 3.0*(fx - fp)/(stp - stx) + dx + dp;
    s = fabs(theta);
    if (s < (t = fabs(dx))) s = t;
    if (s < (t = fabs(dp))) s = t;
    /* The case GAMMA = 0 only arises if the cubic does not tend to
       infinity in the direction of the step. */
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
      /* A minimizer has been bracketed.  If the cubic step is closer to
	 STP than the secant step, the cubic step is taken, otherwise the
	 secant step is taken. */
      stpf = fabs(stpc - stp) < fabs(stpq - stp) ? stpc : stpq;
      t = stp + 0.66*(sty - stp);
      if (stp > stx ? stpf > t : stpf < t) stpf = t;
    } else {
      /* A minimizer has not been bracketed. If the cubic step is farther
	 from stp than the secant step, the cubic step is taken, otherwise
	 the secant step is taken. */
      stpf = fabs(stpc - stp) > fabs(stpq - stp) ? stpc : stpq;
      if (stpf > stpmax) stpf = stpmax;
      if (stpf < stpmin) stpf = stpmin;
    }
  } else {
    /* Fourth case.  A lower function value, derivatives of the same sign,
       and the magnitude of the derivative does not decrease.  If the
       minimum is not bracketed, the step is either STPMIN or STPMAX,
       otherwise the cubic step is taken. */
    info = 4;
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

  /* Update the interval which contains a minimizer and guess for next
     step. */
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
  return info;
}

/*---------------------------------------------------------------------------*/
