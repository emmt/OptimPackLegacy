/*
 * op_utils --
 *
 * Utilities routines for OptimPack library.
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (c) 2003, 2016 Éric Thiébaut.
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
#include "optimpack.h"

/*---------------------------------------------------------------------------*/

int op_error(char *buf, const char *errmsg)
{
  if (buf) {
    if (! errmsg) errmsg = "unknown error";
    strncpy(buf, errmsg, OP_MSG_LEN);
    buf[OP_MSG_LEN] = 0;
  }
  return OP_ERROR;
}

void op_mcopy(const char *msg, char *buf)
{
  if (buf) {
    if (msg) {
      strncpy(buf, msg, OP_MSG_LEN);
      buf[OP_MSG_LEN] = 0;
    } else {
      buf[0] = 0;
    }
  }
}

/*---------------------------------------------------------------------------*/
/* APPLY BOUND CONSTRAINTS */

#define ACTIVE_LO_GRD(x, lo, g) ((x) > (lo) || (g) < zero)
#define ACTIVE_HI_GRD(x, hi, g) ((x) < (hi) || (g) > zero)
#define ACTIVE_LO_DIR(x, lo, d) ((x) > (lo) || (d) > zero)
#define ACTIVE_HI_DIR(x, hi, d) ((x) < (hi) || (d) < zero)

op_integer_t op_bounds_check(op_integer_t n, const double xmin[],
			     const double xmax[])
{
  op_integer_t i;
  if (xmin && xmax) {
    for (i=0 ; i<n ; ++i) {
      if (xmin[i] > xmax[i]) return i;
    }
  }
  return -1;
}

void op_bounds_apply(op_integer_t n, double x[],
		     const double xmin[], const double xmax[])
{
  op_integer_t i;
  if (xmin) {
    if (xmax) {
      /* Lower _and_ upper bounds. */
      for (i=0 ; i<n ; ++i) {
	/**/ if (x[i] < xmin[i]) x[i] = xmin[i];
	else if (x[i] > xmax[i]) x[i] = xmax[i];
      }
    } else {
      /* Only lower bounds. */
      for (i=0 ; i<n ; ++i) {
	if (x[i] < xmin[i]) x[i] = xmin[i];
      }
    }
  } else if (xmax) {
    /* Only upper bounds. */
    for (i=0 ; i<n ; ++i) {
      if (x[i] > xmax[i]) x[i] = xmax[i];
    }
  }
}

void op_bounds_active(op_integer_t n, op_logical_t active[],
		      const double x[], const double g[],
		      const double xmin[], const double xmax[])
{
  const double zero = 0.0;
  op_integer_t i;

  if (xmin) {
    if (xmax) {
      /* Lower _and_ upper bounds. */
      for (i=0 ; i<n ; ++i) {
	active[i] = (ACTIVE_LO_GRD(x[i], xmin[i], g[i]) &&
		     ACTIVE_HI_GRD(x[i], xmax[i], g[i]));
      }
    } else {
      /* Only lower bounds. */
      for (i=0 ; i<n ; ++i) {
	active[i] = ACTIVE_LO_GRD(x[i], xmin[i], g[i]);
      }
    }
  } else if (xmax) {
    /* Only upper bounds. */
    for (i=0 ; i<n ; ++i) {
      active[i] = ACTIVE_HI_GRD(x[i], xmax[i], g[i]);
    }
  }
}

void op_lower_bound_apply(op_integer_t n, double x[], double xmin)
{
  op_integer_t i;
  for (i=0 ; i<n ; ++i) if (x[i] < xmin) x[i] = xmin;
}

void op_lower_bound_active(op_integer_t n, op_logical_t active[],
			   const double x[], const double g[],
			   double xmin)
{
  const double zero = 0.0;
  op_integer_t i;
  for (i=0 ; i<n ; ++i) active[i] = ACTIVE_LO_GRD(x[i], xmin, g[i]);
}

void op_upper_bound_apply(op_integer_t n, double x[], double xmax)
{
  op_integer_t i;
  for (i=0 ; i<n ; ++i) if (x[i] > xmax) x[i] = xmax;
}

void op_upper_bound_active(op_integer_t n, op_logical_t active[],
			   const double x[], const double g[],
			   double xmax)
{
  const double zero = 0.0;
  op_integer_t i;
  for (i=0 ; i<n ; ++i) active[i] = ACTIVE_HI_GRD(x[i], xmax, g[i]);
}

void op_interval_apply(op_integer_t n, double x[], double a, double b)
{
  op_integer_t i;
  if (a > b) { double c=a; a=b; b=c; }
  for (i=0 ; i<n ; ++i) {
    /**/ if (x[i] < a) x[i] = a;
    else if (x[i] > b) x[i] = b;
  }
}

void op_interval_active(op_integer_t n, op_logical_t active[],
			const double x[], const double g[],
			double a, double b)
{
  const double zero = 0.0;
  op_integer_t i;
  if (a > b) { double c=a; a=b; b=c; }
  for (i=0 ; i<n ; ++i) {
    active[i] = (ACTIVE_LO_GRD(x[i], a, g[i]) && ACTIVE_HI_GRD(x[i], b, g[i]));
  }
}

/*---------------------------------------------------------------------------*/
/* VARIANTS OF LEVEL 1 BLAS ROUTINES */

double op_dnrm2(op_integer_t n, const double x[])
{
  const double one = 1.0, zero = 0.0;
  if (n > 1) {
    op_integer_t i;
    double ssq = zero, scale = zero;
    for (i=0 ; i<n ; ++i) {
      if (x[i]) {
	double absxi = fabs(x[i]);
	if (scale < absxi) {
	  double tmp = scale/absxi;
	  ssq = one + ssq*tmp*tmp;
	  scale = absxi;
	} else {
	  double tmp = absxi/scale;
	  ssq += tmp*tmp;
	}
      }
    }
    return scale*sqrt(ssq);
  } else if (n == 1) {
    return fabs(x[0]);
  }
  return zero;
}

void op_dscal(op_integer_t n, double a, double x[])
{
  op_integer_t i;
  if (a == 0.0) {
    memset(x, 0, n*sizeof(double));
  } else if (a == -1.0) {
    for (i=0 ; i<n ; ++i) x[i] -= x[i];
  } else if (a != 1.0) {
    for (i=0 ; i<n ; ++i) x[i] *= a;
  }
}

void op_dcopy(op_integer_t n, const double x[], double y[])
{
  memcpy(y, x, n*sizeof(double));
}

void op_dcopy_active(op_integer_t n, const double x[], double y[],
		     const op_logical_t active[])
{
  if (active) {
    const double zero = 0.0;
    op_integer_t i;
    for (i=0 ; i<n ; ++i) y[i] = (active[i] ? x[i] : zero);
  } else {
    memcpy(y, x, n*sizeof(double));
  }
}

void op_daxpy(op_integer_t n, double a, const double x[], double y[])
{
  op_integer_t i;
  if (a == 1.0) {
    for (i=0 ; i<n ; ++i) y[i] += x[i];
  } else if (a == -1.0) {
    for (i=0 ; i<n ; ++i) y[i] -= x[i];
  } else if (a != 0.0) {
    for (i=0 ; i<n ; ++i) y[i] += a*x[i];
  }
}

void op_daxpy_active(op_integer_t n, double a, const double x[],
		     double y[], const op_logical_t active[])
{
  op_integer_t i;
  if (active) {
    if (a == 1.0) {
      for (i=0 ; i<n ; ++i) if (active[i]) y[i] += x[i];
    } else if (a == -1.0) {
      for (i=0 ; i<n ; ++i) if (active[i]) y[i] -= x[i];
    } else if (a != 0.0) {
      for (i=0 ; i<n ; ++i) if (active[i]) y[i] += a*x[i];
    }
  } else {
    if (a == 1.0) {
      for (i=0 ; i<n ; ++i) y[i] += x[i];
    } else if (a == -1.0) {
      for (i=0 ; i<n ; ++i) y[i] -= x[i];
    } else if (a != 0.0) {
      for (i=0 ; i<n ; ++i) y[i] += a*x[i];
    }
  }
}

double op_ddot(op_integer_t n, const double x[], const double y[])
{
  op_integer_t i;
  double s = 0.0;
  for (i=0 ; i<n ; ++i) s += x[i]*y[i];
  return s;
}

double op_ddot_active(op_integer_t n, const double x[], const double y[],
		      const op_logical_t active[])
{
  op_integer_t i;
  double s = 0.0;
  if (active) {
    for (i=0 ; i<n ; ++i) if (active[i]) s += x[i]*y[i];
  } else {
    for (i=0 ; i<n ; ++i) s += x[i]*y[i];
  }
  return s;
}

/*---------------------------------------------------------------------------*/
/* YORICK-LIKE ROUTINES */

int op_anyof(op_integer_t n, const double x[])
{
  op_integer_t i;
  for (i=0 ; i<n ; ++i) if (x[i]) return 1;
  return 0;
}

int op_noneof(op_integer_t n, const double x[])
{
  op_integer_t i;
  for (i=0 ; i<n ; ++i) if (x[i]) return 0;
  return 1;
}

int op_allof(op_integer_t n, const double x[])
{
  op_integer_t i;
  for (i=0 ; i<n ; ++i) if (! x[i]) return 0;
  return 1;
}

/*---------------------------------------------------------------------------*/
