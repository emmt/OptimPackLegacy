/*
 * opl_yorick.c --
 *
 * Yorick interface for OptimPackLegacy library.
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
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

/* Yorick API. */
#include <pstdlib.h>
#include <play.h>
#include <yapi.h>

#include "optimpacklegacy.h"

/* Define some macros to get rid of some GNU extensions when not compiling
   with GCC. */
#if ! (defined(__GNUC__) && __GNUC__ > 1)
#   define __attribute__(x)
#   define __inline__
#   define __FUNCTION__        ""
#   define __PRETTY_FUNCTION__ ""
#endif

#define TRUE  1
#define FALSE 0

PLUG_API void y_error(const char *) __attribute__ ((noreturn));
static void error(const char *, ...) __attribute__ ((noreturn));

/*---------------------------------------------------------------------------*/
/* PRIVATE FUNCTIONS AND DATA */

/* Indexes for fast parsing of keywords/members. */
static long index_of_dims = -1;
static long index_of_size = -1;
static long index_of_mem = -1;
static long index_of_task = -1;
static long index_of_evaluations = -1;
static long index_of_iterations = -1;
static long index_of_restarts = -1;
static long index_of_step = -1;
static long index_of_gnorm = -1;
static long index_of_fmin = -1;
static long index_of_fatol = -1;
static long index_of_frtol = -1;
static long index_of_sftol = -1;
static long index_of_sgtol = -1;
static long index_of_sxtol = -1;
static long index_of_delta = -1;
static long index_of_epsilon = -1;
static long index_of_status = -1;
static long index_of_reason = -1;

static void
error(const char* format, ...)
{
  char buffer[256];
  const size_t size = sizeof(buffer);
  va_list ap;
  size_t nchars;

  va_start(ap, format);
  nchars = vsnprintf(buffer, size, format, ap);
  va_end(ap);
  if (nchars >= size) {
    buffer[size - 6] = '[';
    buffer[size - 5] = '.';
    buffer[size - 4] = '.';
    buffer[size - 3] = '.';
    buffer[size - 2] = ']';
    buffer[size - 1] = '\0';
  }
  y_error(buffer);
}

/* Push a string on top of the stack. */
static void
push_string(const char* str)
{
  ypush_q(NULL)[0] = p_strcpy(str);
}

/* Define a Yorick global symbol with a long value. */
static void
define_long_const(const char* name, long value)
{
  ypush_long(value);
  yput_global(yget_global(name, 0), 0);
  yarg_drop(1);
}

/* Get Yorick type name. */
static const char*
type_name(int type)
{
  switch (type) {
  case Y_CHAR: return "char";
  case Y_SHORT: return "short";
  case Y_INT: return "int";
  case Y_LONG: return "long";
  case Y_FLOAT: return "float";
  case Y_DOUBLE: return "double";
  case Y_COMPLEX: return "complex";
  case Y_STRING: return "string";
  case Y_POINTER: return "pointer";
  case Y_STRUCT: return "struct";
  case Y_RANGE: return "range";
  case Y_VOID: return "void";
  case Y_FUNCTION: return "function";
  case Y_BUILTIN: return "builtin";
  case Y_STRUCTDEF: return "structdef";
  case Y_STREAM: return "stream";
  case Y_OPAQUE: return "opaque";
  default: return "<unknown_type_number>";
  }
}

/* Get an array of given type and dimensions. */
static void*
get_array(int iarg, int type, const long dims[],
          const char* name, int nil_ok)
{
  int argtype = yarg_typeid(iarg);
  if (argtype == type) {
    long argdims[Y_DIMSIZE];
    long i, rank;
    void* ptr = ygeta_any(iarg, NULL, argdims, NULL);
    if (argdims[0] != (rank = dims[0])) {
      error("bad number of dimensions for argument `%s`", name);
    }
    for (i = 1; i <= rank; ++i) {
      if (argdims[i] != dims[i]) {
        error("bad dimension(s) for argument `%s`", name);
      }
    }
    return ptr;
  } else if (nil_ok && argtype == Y_VOID) {
    return NULL;
  } else {
    error("bad data type for argument `%s`", name);
  }
  error("argument `%s` must be a %ld-D array of `%s`%s",
        name, dims[0], type_name(type), (nil_ok ? " our nil" : ""));
  return NULL;
}

/* Get a dimension list from the stack. */
static long
get_dims(int iarg, long dims[])
{
  int type = yarg_typeid(iarg);
  if (type <= Y_LONG) {
    long j, n, ntot,  rank;
    long* arr = ygeta_l(iarg, &n, dims);
    if (dims[0] == 0 || (dims[0] == 1 && arr[0] == n - 1)) {
      if (dims[0] == 0) {
        /* Got a scalar dimension. */
        rank = 1;
        arr -= 1; /* offset for indexing with j */
      } else {
        /* Got a vector of dimensions. */
        if (n > Y_DIMSIZE) {
          y_error("too many dimensions");
        }
        rank = n - 1;
      }
      ntot = 1;
      dims[0] = rank;
      for (j = 1; j <= rank; ++j) {
        long len = arr[j];
        if (len < 1) {
          y_error("invalid dimension(s)");
        }
        dims[j] = len;
        ntot *= len;
      }
      return ntot;
    }
  } else if (type == Y_VOID) {
    dims[0] = 0;
    return 1;
  }
  y_error("invalid dimension list");
  return -1;
}

/* Function to call before critical code to lower the risk of being
   interrupted. */
static void
critical()
{
  if (p_signalling) {
    p_abort();
  }
}


/*---------------------------------------------------------------------------*/
/* IMPLEMENTATION OF VMLMB INSTANCE */

/* Virtual methods for a VMLMB object. */
static void vmlmb_free(void* ptr);
static void vmlmb_print(void* ptr);
static void vmlmb_eval(void* ptr, int argc);
static void vmlmb_extract(void* ptr, char* name);

/* VMLMB handle type. */
static struct y_userobj_t vmlmb_type = {
  "VMLMB workspace", vmlmb_free, vmlmb_print, vmlmb_eval, vmlmb_extract, NULL
};

/* VMLMB handle instance. (FIXME: save the dimensions as well). */
typedef struct _vmlmb_object {
  opl_vmlmb_workspace_t* ws;
  long n, m;
  long dims[Y_DIMSIZE];
} vmlmb_object_t;

static void
vmlmb_free(void* ptr)
{
  opl_vmlmb_destroy(((vmlmb_object_t*)ptr)->ws);
}

static void
vmlmb_print(void* ptr)
{
  char buffer[128];
  vmlmb_object_t* obj = (vmlmb_object_t*)ptr;
  long i, rank = obj->dims[0];
  sprintf(buffer, "%s with size=%ld, dims=[",
          vmlmb_type.type_name, obj->n);
  y_print(buffer, FALSE);
  for (i = 0; i <= rank; ++i) {
    sprintf(buffer, (i == 0 ? "%ld" : ",%ld"), obj->dims[i]);
    y_print(buffer, FALSE);
  }
  sprintf(buffer, "], mem=%ld, task=%d",
          obj->m, opl_vmlmb_get_task(obj->ws));
  y_print(buffer, TRUE);
}

static void
vmlmb_eval(void* ptr, int argc)
{
  ypush_nil();
}

static void
vmlmb_extract(void* ptr, char* name)
{
  vmlmb_object_t* obj = (vmlmb_object_t*)ptr;
  long index = yget_global(name, 0);
  if (index == index_of_iterations) {
    ypush_long(opl_vmlmb_get_iterations(obj->ws));
  } else if (index == index_of_evaluations) {
    ypush_long(opl_vmlmb_get_evaluations(obj->ws));
  } else if (index == index_of_restarts) {
    ypush_long(opl_vmlmb_get_restarts(obj->ws));
  } else if (index == index_of_task) {
    ypush_long(opl_vmlmb_get_task(obj->ws));
  } else if (index == index_of_step) {
    ypush_double(opl_vmlmb_get_step(obj->ws));
  } else if (index == index_of_gnorm) {
    ypush_double(opl_vmlmb_get_gnorm(obj->ws));
  } else if (index == index_of_fmin) {
    ypush_double(opl_vmlmb_get_fmin(obj->ws));
  } else if (index == index_of_fatol) {
    ypush_double(opl_vmlmb_get_fatol(obj->ws));
  } else if (index == index_of_frtol) {
    ypush_double(opl_vmlmb_get_frtol(obj->ws));
  } else if (index == index_of_sftol) {
    ypush_double(opl_vmlmb_get_sftol(obj->ws));
  } else if (index == index_of_sgtol) {
    ypush_double(opl_vmlmb_get_sgtol(obj->ws));
  } else if (index == index_of_sxtol) {
    ypush_double(opl_vmlmb_get_sxtol(obj->ws));
  } else if (index == index_of_delta) {
    ypush_double(opl_vmlmb_get_delta(obj->ws));
  } else if (index == index_of_epsilon) {
    ypush_double(opl_vmlmb_get_epsilon(obj->ws));
  } else if (index == index_of_size) {
    ypush_long(obj->n);
  } else if (index == index_of_mem) {
    ypush_long(obj->m);
  } else if (index == index_of_dims) {
    long i, rank = obj->dims[0];
    long dims[2];
    long* arr;
    dims[0] = 1;
    dims[1] = rank + 1;
    arr = ypush_l(dims);
    for (i = 0; i <= rank; ++i) {
      arr[i] = obj->dims[i];
    }
  } else if (index == index_of_status) {
    ypush_long(opl_vmlmb_get_status(obj->ws));
  } else if (index == index_of_reason) {
    push_string(opl_vmlmb_get_reason(obj->ws));
  } else {
    y_error("unknown member");
  }
}

static vmlmb_object_t*
get_vmlmb(int iarg)
{
  return (vmlmb_object_t*)yget_obj(iarg, &vmlmb_type);
}

/*---------------------------------------------------------------------------*/
/* BUILTIN FUNCTIONS */

void
Y_opl_vmlmb_create(int argc)
{
  long i, rank, m = -1, n = -1;
  long dims[Y_DIMSIZE];
  vmlmb_object_t* obj;
  int iarg;
  int fmin_iarg = -1;
  int fatol_iarg = -1;
  int frtol_iarg = -1;
  int sftol_iarg = -1;
  int sgtol_iarg = -1;
  int sxtol_iarg = -1;
  int delta_iarg = -1;
  int epsilon_iarg = -1;

  /* Parse arguments. */
  for (iarg = argc - 1; iarg >= 0; --iarg) {
    long index = yarg_key(iarg);
    if (index < 0) {
      /* Positional argument. */
      if (n == -1) {
        n = get_dims(iarg, dims);
      } else if (m == -1) {
        m = ygets_l(iarg);
        if (m <= 0) {
          y_error("invalid number of steps to memorize");
        }
        if (m > n) {
          m = n;
        }
      } else {
        y_error("too many arguments");
      }
    } else {
      /* Keyword argument (skip its value). */
      --iarg;
      if (index == index_of_fmin) {
        fmin_iarg = iarg;
      } else if (index == index_of_fatol) {
        fatol_iarg = iarg;
      } else if (index == index_of_frtol) {
        frtol_iarg = iarg;
      } else if (index == index_of_sftol) {
        sftol_iarg = iarg;
      } else if (index == index_of_sgtol) {
        sgtol_iarg = iarg;
      } else if (index == index_of_sxtol) {
        sxtol_iarg = iarg;
      } else if (index == index_of_delta) {
        delta_iarg = iarg;
      } else if (index == index_of_epsilon) {
        epsilon_iarg = iarg;
      } else {
        y_error("unsupported keyword");
      }
    }
  }
  if (n == -1) {
    y_error("missing dimension list of variables");
  }
  if (m == -1) {
    y_error("missing number of steps to memorize");
  }

  /* Create VMLMB instance. */
  obj = (vmlmb_object_t*)ypush_obj(&vmlmb_type, sizeof(vmlmb_object_t));
  critical();
  obj->ws = opl_vmlmb_create(n, m);
  if (obj->ws == NULL) {
    if (errno == ENOMEM) {
      y_error("insufficient memory");
    } else {
      y_error("unknown error");
    }
  }
  obj->n = n;
  obj->m = m;
  rank = dims[0];
  for (i = 0; i <= rank; ++i) {
    obj->dims[i] = dims[i];
  }

  /* Configure VMLMB instance (adding +1 to iarg's because an element has been
     pushed on top of the stack).. */
# define SET_ATTRIBUTE(name, invalid)                           \
  if (name##_iarg >= 0 && ! yarg_nil(name##_iarg + 1)) {        \
    double value = ygets_d(name##_iarg + 1);                    \
    if ((invalid) ||                                            \
        opl_vmlmb_set_##name(obj->ws, value) != OPL_SUCCESS) {  \
      y_error("invalid value for `" #name "`");                 \
    }                                                           \
  }
  SET_ATTRIBUTE(fmin, FALSE);
  SET_ATTRIBUTE(fatol, value < 0);
  SET_ATTRIBUTE(frtol, value < 0);
  SET_ATTRIBUTE(sftol, value <= 0 || value >= 1);
  SET_ATTRIBUTE(sgtol, value <= 0 || value >= 1);
  SET_ATTRIBUTE(sxtol, value <= 0 || value >= 1);
  SET_ATTRIBUTE(delta, value < 0);
  SET_ATTRIBUTE(epsilon, value < 0);
# undef SET_ATTRIBUTE

}

void
Y_opl_vmlmb_configure(int argc)
{
  vmlmb_object_t* obj = NULL;
  int iarg;
  int ndrop = 0;
  int fmin_iarg = -1;
  int fatol_iarg = -1;
  int frtol_iarg = -1;
  int sftol_iarg = -1;
  int sgtol_iarg = -1;
  int sxtol_iarg = -1;
  int delta_iarg = -1;
  int epsilon_iarg = -1;

  /* Parse arguments. */
  for (iarg = argc - 1; iarg >= 0; --iarg) {
    long index = yarg_key(iarg);
    if (index < 0) {
      /* Positional argument. */
      ndrop += 1;
      if (obj == NULL) {
        obj = get_vmlmb(iarg);
        ndrop = 0;
      } else {
        y_error("too many arguments");
      }
    } else {
      /* Keyword argument (skip its value). */
      ndrop += 2;
      --iarg;
      if (index == index_of_fmin) {
        fmin_iarg = iarg;
      } else if (index == index_of_fatol) {
        fatol_iarg = iarg;
      } else if (index == index_of_frtol) {
        frtol_iarg = iarg;
      } else if (index == index_of_sftol) {
        sftol_iarg = iarg;
      } else if (index == index_of_sgtol) {
        sgtol_iarg = iarg;
      } else if (index == index_of_sxtol) {
        sxtol_iarg = iarg;
      } else if (index == index_of_delta) {
        delta_iarg = iarg;
      } else if (index == index_of_epsilon) {
        epsilon_iarg = iarg;
      } else {
        y_error("unsupported keyword");
      }
    }
  }
  if (obj == NULL) {
    y_error("missing VMLMB workspace");
  }

  /* Set attributes. */
# define SET_ATTRIBUTE(name, invalid)                           \
  if (name##_iarg >= 0 && ! yarg_nil(name##_iarg)) {            \
    double value = ygets_d(name##_iarg);                        \
    if ((invalid) ||                                            \
        opl_vmlmb_set_##name(obj->ws, value) != OPL_SUCCESS) {  \
      y_error("invalid value for `" #name "`");                 \
    }                                                           \
  }
  SET_ATTRIBUTE(fmin, FALSE);
  SET_ATTRIBUTE(fatol, value < 0);
  SET_ATTRIBUTE(frtol, value < 0);
  SET_ATTRIBUTE(sftol, value <= 0 || value >= 1);
  SET_ATTRIBUTE(sgtol, value <= 0 || value >= 1);
  SET_ATTRIBUTE(sxtol, value <= 0 || value >= 1);
  SET_ATTRIBUTE(delta, value < 0);
  SET_ATTRIBUTE(epsilon, value < 0);
# undef SET_ATTRIBUTE

  /* Manage to left WS on top of the stack. */
  if (ndrop > 0) {
    yarg_drop(ndrop);
  }
}

void
Y_opl_vmlmb_iterate(int argc)
{
  double f;
  vmlmb_object_t* obj;
  int* isfree;
  double* x;
  double* g;
  double* h;
  opl_task_t task;
  long fref;
  int iarg;

  if (argc < 4 || argc > 6) {
    y_error("expecting between 4 and 6 arguments");
  }
  iarg = argc;
  obj = get_vmlmb(--iarg);
  x = (double*)get_array(--iarg, Y_DOUBLE, obj->dims, "x", FALSE);
  fref = yget_ref(--iarg);
  if (fref < 0) {
    y_error("expecting a simple variable reference for argument `f`");
  }
  f = ygets_d(iarg);
  g = (double*)get_array(--iarg, Y_DOUBLE, obj->dims, "g", FALSE);
  if (argc >= 5) {
    isfree = (int*)get_array(--iarg, Y_INT, obj->dims, "isfree", TRUE);
  } else {
    isfree = NULL;
  }
  if (argc >= 6) {
    h = (double*)get_array(--iarg, Y_DOUBLE, obj->dims, "h", TRUE);
  } else {
    h = NULL;
  }
  task = opl_vmlmb_iterate(obj->ws, x, &f, g, isfree, h);
  ypush_double(f);
  yput_global(fref, 0);
  ypush_long(task);
}

void
Y_opl_vmlmb_restart(int argc)
{
  vmlmb_object_t* obj;

  if (argc != 1) {
    y_error("expecting exactly one argument");
  }
  obj = get_vmlmb(0);
  opl_vmlmb_restart(obj->ws);
  ypush_long(opl_vmlmb_get_task(obj->ws));
}

void
Y_opl_vmlmb_restore(int argc)
{
  double f;
  vmlmb_object_t* obj;
  double* x;
  double* g;
  long fref;
  int iarg;

  if (argc != 4) {
    y_error("expecting exactly 4 arguments");
  }
  iarg = argc;
  obj = get_vmlmb(--iarg);
  x = (double*)get_array(--iarg, Y_DOUBLE, obj->dims, "x", FALSE);
  fref = yget_ref(--iarg);
  if (fref < 0) {
    y_error("expecting a simple variable reference for argument `f`");
  }
  g = (double*)get_array(--iarg, Y_DOUBLE, obj->dims, "g", FALSE);
  opl_vmlmb_restore(obj->ws, x, &f, g);
  ypush_double(f);
  yput_global(fref, 0);
  ypush_long(opl_vmlmb_get_task(obj->ws));
}

void
Y__opl_init(int argc)
{
  /* Define constants. */
#define DEFINE_LONG_CONST(c) define_long_const(#c, c)
  DEFINE_LONG_CONST(OPL_TASK_START);
  DEFINE_LONG_CONST(OPL_TASK_FG);
  DEFINE_LONG_CONST(OPL_TASK_FREEVARS);
  DEFINE_LONG_CONST(OPL_TASK_NEWX);
  DEFINE_LONG_CONST(OPL_TASK_CONV);
  DEFINE_LONG_CONST(OPL_TASK_WARN);
  DEFINE_LONG_CONST(OPL_TASK_ERROR);
#undef DEFINE_LONG_CONST

  /* Define fast keyword/member indexes. */
#define INIT(s) if (index_of_##s == -1L) index_of_##s = yget_global(#s, 0)
  INIT(dims);
  INIT(size);
  INIT(mem);
  INIT(task);
  INIT(evaluations);
  INIT(iterations);
  INIT(restarts);
  INIT(step);
  INIT(gnorm);
  INIT(fmin);
  INIT(fatol);
  INIT(frtol);
  INIT(sftol);
  INIT(sgtol);
  INIT(sxtol);
  INIT(delta);
  INIT(epsilon);
  INIT(status);
  INIT(reason);
#undef INIT

  /* In case of... */
  ypush_nil();
}

