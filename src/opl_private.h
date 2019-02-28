/*
 * opl_private.h --
 *
 * Private definitions for optimization routines implemented in Optimpacklegacy
 * library.
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

#ifndef _OPL_PRIVATE_H
#define _OPL_PRIVATE_H 1

#include "optimpacklegacy.h"

/*---------------------------------------------------------------------------*/
/* USEFUL MACROS */

#ifndef __STRICT_ANSI__
#ifndef _OPL_FORCE_INLINE
#  if defined(_MSC_VER) /* Microsoft Visual Studio */
#    define _OPL_FORCE_INLINE __forceinline
#  elif defined(__GNUC__) /* GNU Compiler */
#    define _OPL_FORCE_INLINE static inline __attribute__((always_inline))
#  elif defined( __cplusplus) /* C++ Compiler */
#    define _OPL_FORCE_INLINE inline
#  endif
#endif
#endif

/* OPL_STRINGIFY takes an argument and wraps it in "" (double quotation
   marks), OPL_CONCAT concatenates two arguments. */
#ifdef __STDC__
# define OPL_STRINGIFY(x)     #x
# define OPL_CONCAT(a,b)      a##b
# define OPL_CONCAT2(a,b)     a##b
# define OPL_CONCAT3(a,b,c)   a##b##c
# define OPL_CONCAT4(a,b,c,d) a##b##c##d
#else
# define OPL_STRINGIFY(x)     "x"
# define OPL_CONCAT(a,b)      a/**/b
# define OPL_CONCAT2(a,b)     a/**/b
# define OPL_CONCAT3(a,b,c)   a/**/b/**/c
# define OPL_CONCAT4(a,b,c,d) a/**/b/**/c/**/d
#endif

/* Computes absolute value: */
#define OPL_ABS(a)   ((a)>=0?(a):-(a))

/* Computes min/max values: */
#define OPL_MIN(a,b) ((a)<=(b)?(a):(b))
#define OPL_MAX(a,b) ((a)>=(b)?(a):(b))

/* Computes minimal number of chunks with M elements
   needed to store N elements: */
#define OPL_HOW_MANY(n, m) (((n)+((m)-1))/(m))

/* Returns N elements rounding up to a multiple of M elements: */
#define OPL_ROUND_UP(n, m) (OPL_HOW_MANY(n, m)*(m))

/* Offset (in bytes) of member M in structure S: */
#define OPL_OFFSET_OF(s, m) ((size_t) &((s *)0)->m)

/* Allocate memory for a given number of elements of a given type. */
#define OPL_NEW(type, number)  (type*)malloc((number)*sizeof(type))

/*---------------------------------------------------------------------------*/
/* STRUCTURES */

/**
 * The context structure.
 */
struct _opl_context {
  const char*  message;     /* Current message. */
  opl_status_t status;      /* Status code. */
  int          syserr;      /* Value of `errno` if status is
                               `OPL_SYSTEM_ERROR` */
  char         buffer[128]; /* Internal buffer */
};

/** Global message for successful operation. */
extern const char* _opl_success_message;

/**
 * Set a context with a permanent message.
 */
#define _OPL_SET_CONTEXT(ctx, stat, mesg) \
  ((ctx)->message = (mesg), (ctx)->status = stat)

#ifdef _OPL_FORCE_INLINE
opl_status_t  _OPL_FORCE_INLINE _opl_success(opl_context_t* ctx)
{
  return _OPL_SET_CONTEXT(ctx, OPL_SUCCESS, _opl_success_message);
}
#endif

/**
 * Workspace data used for Moré & Thuente line search.
 */
struct _opl_csrch_workspace {
  opl_context_t context;
  double        ftol, gtol, xtol;
  double        stpmin, stpmax;
  double        finit;
  double        ginit;
  double        stx, fx, gx;
  double        sty, fy, gy;
  double        stmin, stmax;
  double        width, width1;
  opl_task_t    task;
  int           stage;
  opl_boolean_t brackt;
};

/**
 * The workspace for VMLMB.
 */
struct _opl_vmlmb_workspace {
  opl_csrch_workspace_t lnsrch;
  opl_integer_t n;
  opl_integer_t m;
  opl_integer_t mp;
  opl_integer_t mark;
  opl_integer_t evaluations;
  opl_integer_t iterations;
  opl_integer_t restarts;
  unsigned int flags;
  void (*free)(void*);
  double frtol;
  double fatol;
  double fmin;
  double f0;
  double gd;
  double g0d;
  double stp;
  double delta;
  double epsilon;
  double gnorm;
  double g0norm;
  double* alpha;
  double* rho;
  double* d;
  double** S;
  double** Y;
};

#endif /* _OPL_PRIVATE_H */
