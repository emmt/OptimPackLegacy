/*
 * opl_algebra.c --
 *
 * Basic linear algebra routines for OptimPackLegacy library.
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
#include "opl_private.h"

/*---------------------------------------------------------------------------*/
/* APPLY BOUND CONSTRAINTS */

#define ISFREE_LO_GRD(x, lo, g) ((x) > (lo) || (g) < zero)
#define ISFREE_HI_GRD(x, hi, g) ((x) < (hi) || (g) > zero)
#define ISFREE_LO_DIR(x, lo, d) ((x) > (lo) || (d) > zero)
#define ISFREE_HI_DIR(x, hi, d) ((x) < (hi) || (d) < zero)

opl_integer_t opl_bounds_check(opl_integer_t n, const double xmin[],
                               const double xmax[])
{
  opl_integer_t i;
  if (xmin && xmax) {
    for (i = 0; i < n; ++i) {
      if (xmin[i] > xmax[i]) {
        return i;
      }
    }
  }
  return -1;
}

void opl_bounds_apply(opl_integer_t n, double x[],
                      const double xmin[], const double xmax[])
{
  opl_integer_t i;
  if (xmin) {
    if (xmax) {
      /* Lower _and_ upper bounds. */
      for (i = 0; i < n; ++i) {
	if (x[i] < xmin[i]) {
          x[i] = xmin[i];
        }
	if (x[i] > xmax[i]) {
          x[i] = xmax[i];
        }
      }
    } else {
      /* Only lower bounds. */
      for (i = 0; i < n; ++i) {
	if (x[i] < xmin[i]) {
          x[i] = xmin[i];
        }
      }
    }
  } else if (xmax) {
    /* Only upper bounds. */
    for (i = 0; i < n; ++i) {
      if (x[i] > xmax[i]) {
        x[i] = xmax[i];
      }
    }
  }
}

void opl_bounds_free(opl_integer_t n, opl_logical_t isfree[],
                     const double x[], const double g[],
                     const double xmin[], const double xmax[])
{
  const double zero = 0.0;
  opl_integer_t i;

  if (xmin) {
    if (xmax) {
      /* Lower _and_ upper bounds. */
      for (i = 0; i < n; ++i) {
	isfree[i] = (ISFREE_LO_GRD(x[i], xmin[i], g[i]) &&
		     ISFREE_HI_GRD(x[i], xmax[i], g[i]));
      }
    } else {
      /* Only lower bounds. */
      for (i = 0; i < n; ++i) {
	isfree[i] = ISFREE_LO_GRD(x[i], xmin[i], g[i]);
      }
    }
  } else if (xmax) {
    /* Only upper bounds. */
    for (i = 0; i < n; ++i) {
      isfree[i] = ISFREE_HI_GRD(x[i], xmax[i], g[i]);
    }
  }
}

void opl_lower_bound_apply(opl_integer_t n, double x[], double xmin)
{
  opl_integer_t i;
  for (i = 0; i < n; ++i) {
    if (x[i] < xmin) {
      x[i] = xmin;
    }
  }
}

void opl_lower_bound_free(opl_integer_t n, opl_logical_t isfree[],
                          const double x[], const double g[],
                          double xmin)
{
  const double zero = 0.0;
  opl_integer_t i;
  for (i = 0; i < n; ++i) {
    isfree[i] = ISFREE_LO_GRD(x[i], xmin, g[i]);
  }
}

void opl_upper_bound_apply(opl_integer_t n, double x[], double xmax)
{
  opl_integer_t i;
  for (i = 0; i < n; ++i) {
    if (x[i] > xmax) {
      x[i] = xmax;
    }
  }
}

void opl_upper_bound_free(opl_integer_t n, opl_logical_t isfree[],
                          const double x[], const double g[],
                          double xmax)
{
  const double zero = 0.0;
  opl_integer_t i;
  for (i = 0; i < n; ++i) {
    isfree[i] = ISFREE_HI_GRD(x[i], xmax, g[i]);
  }
}

void opl_interval_apply(opl_integer_t n, double x[], double a, double b)
{
  opl_integer_t i;
  if (a > b) {
    double c = a;
    a = b;
    b = c;
  }
  for (i = 0; i < n; ++i) {
    if (x[i] < a) x[i] = a;
    if (x[i] > b) x[i] = b;
  }
}

void opl_interval_free(opl_integer_t n, opl_logical_t isfree[],
                       const double x[], const double g[],
                       double a, double b)
{
  const double zero = 0.0;
  opl_integer_t i;
  if (a > b) { double c=a; a=b; b=c; }
  for (i = 0; i < n; ++i) {
    isfree[i] = (ISFREE_LO_GRD(x[i], a, g[i]) &&
                 ISFREE_HI_GRD(x[i], b, g[i]));
  }
}

/*---------------------------------------------------------------------------*/
/* VARIANTS OF LEVEL 1 BLAS ROUTINES */

double opl_dnrm2(opl_integer_t n, const double x[])
{
  const double one = 1.0, zero = 0.0;
  if (n > 1) {
    opl_integer_t i;
    double ssq = zero, scale = zero;
    for (i = 0; i < n; ++i) {
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

void opl_dscal(opl_integer_t n, double a, double x[])
{
  opl_integer_t i;
  if (a == 0.0) {
    memset(x, 0, n*sizeof(double));
  } else if (a == -1.0) {
    for (i = 0; i < n; ++i) {
      x[i] -= x[i];
    }
  } else if (a != 1.0) {
    for (i = 0; i < n; ++i) {
      x[i] *= a;
    }
  }
}

void opl_dcopy(opl_integer_t n, const double x[], double y[])
{
  memcpy(y, x, n*sizeof(double));
}

void opl_dcopy_free(opl_integer_t n, const double x[], double y[],
                    const opl_logical_t isfree[])
{
  if (isfree != NULL) {
    const double zero = 0.0;
    opl_integer_t i;
    for (i = 0; i < n; ++i) {
      y[i] = (isfree[i] ? x[i] : zero);
    }
  } else {
    memcpy(y, x, n*sizeof(double));
  }
}

void opl_daxpy(opl_integer_t n, double a, const double x[], double y[])
{
  opl_integer_t i;
  if (a == 1.0) {
    for (i = 0; i < n; ++i) {
      y[i] += x[i];
    }
  } else if (a == -1.0) {
    for (i = 0; i < n; ++i) {
      y[i] -= x[i];
    }
  } else if (a != 0.0) {
    for (i = 0; i < n; ++i) {
      y[i] += a*x[i];
    }
  }
}

void opl_daxpy_free(opl_integer_t n, double a, const double x[],
                    double y[], const opl_logical_t isfree[])
{
  opl_integer_t i;
  if (isfree != NULL) {
    if (a == 1.0) {
      for (i = 0; i < n; ++i) {
        if (isfree[i]) {
          y[i] += x[i];
        }
      }
    } else if (a == -1.0) {
      for (i = 0; i < n; ++i) {
        if (isfree[i]) {
          y[i] -= x[i];
        }
      }
    } else if (a != 0.0) {
      for (i = 0; i < n; ++i) {
        if (isfree[i]) {
          y[i] += a*x[i];
        }
      }
    }
  } else {
    if (a == 1.0) {
      for (i = 0; i < n; ++i) {
        y[i] += x[i];
      }
    } else if (a == -1.0) {
      for (i = 0; i < n; ++i) {
        y[i] -= x[i];
      }
    } else if (a != 0.0) {
      for (i = 0; i < n; ++i) {
        y[i] += a*x[i];
      }
    }
  }
}

double opl_ddot(opl_integer_t n, const double x[], const double y[])
{
  opl_integer_t i;
  double s = 0.0;
  for (i = 0; i < n; ++i) {
    s += x[i]*y[i];
  }
  return s;
}

double opl_ddot_free(opl_integer_t n, const double x[], const double y[],
                     const opl_logical_t isfree[])
{
  opl_integer_t i;
  double s = 0.0;
  if (isfree != NULL) {
    for (i = 0; i < n; ++i) {
      if (isfree[i]) {
        s += x[i]*y[i];
      }
    }
  } else {
    for (i = 0; i < n; ++i) {
      s += x[i]*y[i];
    }
  }
  return s;
}

/*---------------------------------------------------------------------------*/
/* YORICK-LIKE ROUTINES */

int opl_anyof(opl_integer_t n, const double x[])
{
  opl_integer_t i;
  for (i = 0; i < n; ++i) {
    if (x[i]) {
      return 1;
    }
  }
  return 0;
}

int opl_noneof(opl_integer_t n, const double x[])
{
  opl_integer_t i;
  for (i = 0; i < n; ++i) {
    if (x[i]) {
      return 0;
    }
  }
  return 1;
}

int opl_allof(opl_integer_t n, const double x[])
{
  opl_integer_t i;
  for (i = 0; i < n; ++i) {
    if (! x[i]) {
      return 0;
    }
  }
  return 1;
}

/*---------------------------------------------------------------------------*/
