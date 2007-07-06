/*
 * op_idl.c --
 *
 *	Routines to interface OptimPack into IDL.
 *
 *-----------------------------------------------------------------------------
 *
 *      Copyright (c) 2003, Eric Thiébaut.
 *
 *	This file is part of OptimPack.
 *
 *	OptimPack is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU General Public License
 *	version 2 as published by the Free Software Foundation.
 *
 *	OptimPack is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU General Public License for more details.
 *
 *	You should have received a copy of the GNU General Public
 *	License along with OptimPack (file "COPYING" in the top source
 *	directory); if not, write to the Free Software Foundation,
 *	Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
 *
 *-----------------------------------------------------------------------------
 *
 *	$Id: op_idl.c,v 1.1 2007/07/06 23:27:02 eric Exp eric $
 *	$Log: op_idl.c,v $
 *	Revision 1.1  2007/07/06 23:27:02  eric
 *	Initial revision
 *
 *-----------------------------------------------------------------------------
 */

#include "op_idl.h"
#include "optimpack.h"

/*---------------------------------------------------------------------------*/
/* MACROS */

/* Macros to parse IDL SIZE array:
   [NDIMS, DIM_1, ..., DIM_NDIMS, TYPE, NELEMENTS] */
#define GET_TYPE(s)      ((s)[*(s) + 1])
#define GET_NDIMS(s)     ((s)[0])
#define GET_NELEMENTS(s) ((s)[*(s) + 2])
#define IS_SCALAR(s)     (GET_NDIMS(s) == 0)

static const char *last_error = 0;
static const char *get_error(int reset);
static void save_error(const char *msg);
static void reset_error();

/* Fetch scalar argument. */
static idl_long_t   get_integer(void *size_ptr, void *addr);
static idl_double_t get_real(void *size_ptr, void *addr);

/* Fetch array argument. */
static idl_byte_t     *get_byte_array(void *size, void *addr,
				      idl_long_t *nelements);
static idl_int_t      *get_int_array(void *size, void *addr,
				     idl_long_t *nelements);
static idl_long_t     *get_long_array(void *size, void *addr,
				      idl_long_t *nelements);
static idl_long64_t   *get_long64_array(void *size, void *addr,
					idl_long_t *nelements);
static idl_float_t    *get_float_array(void *size, void *addr,
				       idl_long_t *nelements);
static idl_double_t   *get_double_array(void *size, void *addr,
					idl_long_t *nelements);
static idl_complex_t  *get_complex_array(void *size, void *addr,
					 idl_long_t *nelements);
static idl_dcomplex_t *get_dcomplex_array(void *size, void *addr,
					  idl_long_t *nelements);

/*---------------------------------------------------------------------------*/

const char *op_idl_last_error(int argc, void *argv[])
{
  return (last_error ? last_error : "no error");
}

/* Macros to unpack workspace arrays for VMLMB. */
#define HOW_MANY(a,b)           OP_HOW_MANY(a,b)
#define ROUND_UP(a,b)           OP_ROUND_UP(a,b)

#define VMLMB_CSAVE_NUMBER      OP_VMLMB_CSAVE_NUMBER
#define VMLMB_ISAVE_NUMBER      OP_VMLMB_ISAVE_NUMBER
#define VMLMB_DSAVE_NUMBER(n,m) OP_VMLMB_DSAVE_NUMBER(n,m)
#define VMLMB_CSAVE_OFFSET 0
#define VMLMB_ISAVE_OFFSET ROUND_UP(VMLMB_CSAVE_OFFSET + 		   \
				    VMLMB_CSAVE_NUMBER*sizeof(char),	   \
				    sizeof(idl_long_t))
#define VMLMB_DSAVE_OFFSET ROUND_UP(VMLMB_ISAVE_OFFSET +		   \
				    VMLMB_ISAVE_NUMBER*sizeof(idl_long_t), \
				    sizeof(idl_double_t))
#define VMLMB_CSAVE(ws) ((char         *)((char *)(ws) + VMLMB_CSAVE_OFFSET))
#define VMLMB_ISAVE(ws) ((idl_long_t   *)((char *)(ws) + VMLMB_ISAVE_OFFSET))
#define VMLMB_DSAVE(ws) ((idl_double_t *)((char *)(ws) + VMLMB_DSAVE_OFFSET))
#define VMLMB_WS_NUMBER(n,m) (VMLMB_DSAVE_NUMBER(n,m) +			    \
			      HOW_MANY(VMLMB_DSAVE_OFFSET, sizeof(double)))

idl_long_t op_idl_vmlmb_first(int argc, void *argv[])
{
  idl_long_t   n, m, task, nws;
  idl_double_t fmin, frtol, fatol, sftol, sgtol, sxtol;
  idl_long_t *isave;
  idl_double_t *dsave;
  void *ws;
  char *csave;

  /* Parse arguments (N, M, WS). */
  if (argc != 6) {
    last_error = "expecting exactly 6 arguments";
    return -1;
  }
  last_error = 0;
  n = get_integer(argv[0], argv[1]);
  m = get_integer(argv[2], argv[3]);
  ws = get_double_array(argv[4], argv[5], &nws);
  if (last_error) return -1;

  /* Check worspace size and unpack workspace arrays and real parameters. */
  if (nws < 6 /* because of packed real parameters */ ||
      nws < VMLMB_WS_NUMBER(n, m)) {
    if (! last_error) last_error = "insufficient space in workspace array";
    return -1;
  }
  fmin  = ((double *)ws)[0];
  frtol = ((double *)ws)[1];
  fatol = ((double *)ws)[2];
  sftol = ((double *)ws)[3];
  sgtol = ((double *)ws)[4];
  sxtol = ((double *)ws)[5];

#if 0
  fprintf(stderr, "NSW   = %d\n", (int)nws);
  fprintf(stderr, "N     = %d\n", (int)n);
  fprintf(stderr, "M     = %d\n", (int)m);
  fprintf(stderr, "FMIN  = %g\n", fmin);
  fprintf(stderr, "FRTOL = %g\n", frtol);
  fprintf(stderr, "FATOL = %g\n", fatol);
  fprintf(stderr, "SFTOL = %g\n", sftol);
  fprintf(stderr, "SGTOL = %g\n", sgtol);
  fprintf(stderr, "SXTOL = %g\n", sxtol);
#endif

  csave = VMLMB_CSAVE(ws);
  isave = VMLMB_ISAVE(ws);
  dsave = VMLMB_DSAVE(ws);

  /* Call op_vmlmb_first and check returned value. */
  task = op_vmlmb_first(n, m,
			fmin, fatol, frtol, sftol, sgtol, sxtol,
			csave, isave, dsave);
  if (task != OP_TASK_FG) {
    save_error(csave);
    return -1;
  }
  return task;
}

idl_long_t op_idl_vmlmb_next(int argc, void *argv[])
{
  idl_long_t *size, task;
  idl_double_t *x, *f, *g, *h, *dsave;
  idl_byte_t *active;
  idl_long_t n, m, number, *isave;
  char *csave;
  void *addr;

  /* Parse arguments (X, F, G, WS[, ACTIVE [, H]]). */
  if (argc != 8 && argc != 10 && argc != 12) {
    last_error = "expecting 8, 10 or 12 arguments";
    return -1;
  }

  /* Get workspace and make sure it is large enough to store parameters N
     and M and, then large enough to store all workspace arrays. */
  size = (idl_long_t *)argv[6];
  addr = argv[7];
  if (GET_TYPE(size) != IDL_TYP_DOUBLE) {
    last_error = "workspace WS should be an array of double's";
    return -1;
  }
  number = GET_NELEMENTS(size);
  if (number*sizeof(idl_double_t) < VMLMB_DSAVE_OFFSET) {
  bad_workspace_size:
    last_error = "insufficient space in workspace array";
    return -1;
  }
  csave = VMLMB_CSAVE(addr);
  isave = VMLMB_ISAVE(addr);
  dsave = VMLMB_DSAVE(addr);
  if ((m = isave[4]) <= 0 || (n = isave[5]) <= 0) {
    last_error = "corrupted workspace array";
    return -1;
  }
  if (number < VMLMB_WS_NUMBER(n, m)) goto bad_workspace_size;

  /* Get parameters X. */
  size = (idl_long_t *)argv[0];
  addr = argv[1];
  if (GET_TYPE(size) != IDL_TYP_DOUBLE) {
    last_error = "parameters X should be an array of DOUBLE's";
    return -1;
  }
  if (GET_NELEMENTS(size) != n) {
    last_error = "bad number of elements for parameter array X";
    return -1;
  }
  x = (idl_double_t *)addr;

  /* Get function value F. */
  size = (idl_long_t *)argv[2];
  addr = argv[3];
  if (GET_NDIMS(size) || GET_TYPE(size) != IDL_TYP_DOUBLE) {
    last_error = "function value F should be a DOUBLE scalar";
    return -1;
  }
  f = (idl_double_t *)addr;

  /* Get gradient G. */
  size = (idl_long_t *)argv[4];
  addr = argv[5];
  if (GET_TYPE(size) != IDL_TYP_DOUBLE) {
    last_error = "gradient G should be an array of DOUBLE's";
    return -1;
  }
  if (GET_NELEMENTS(size) != n) {
    last_error = "bad number of elements for gradient array G";
    return -1;
  }
  g = (idl_double_t *)addr;

  /* Get optional ACTIVE array. */
  if (argc >= 10) {
    size = (idl_long_t *)argv[8];
    addr = argv[9];
    switch (GET_TYPE(size)) {
    case IDL_TYP_BYTE:
      if (GET_NELEMENTS(size) != n) {
	last_error = "bad number of elements for array ACTIVE";
	return -1;
      }
      active = (idl_byte_t *)addr;
      break;
    case IDL_TYP_UNDEF:
      active = (idl_byte_t *)0;
      break;
    default:
      last_error = "ACTIVE should be an array of BYTE's or undefined";
      return -1;
    }
  } else {
    active = (idl_byte_t *)0;
  }

  /* Get optional H array. */
  if (argc >= 12) {
    size = (idl_long_t *)argv[10];
    addr = argv[11];
    switch (GET_TYPE(size)) {
    case IDL_TYP_DOUBLE:
      if (GET_NELEMENTS(size) != n) {
	last_error = "bad number of elements for array H";
	return -1;
      }
      h = (idl_double_t *)addr;
      break;
    case IDL_TYP_UNDEF:
      h = (idl_double_t *)0;
      break;
    default:
      last_error = "H should be an array of DOUBLE's or undefined";
      return -1;
    }
  } else {
    h = (idl_double_t *)0;
  }

  /* Call op_vmlmb_next and check returned value. */
  task = op_vmlmb_next(x, f, g, active, h, csave, isave, dsave);
  if (task == OP_TASK_ERROR) {
    save_error(csave);
  } else {
    last_error = 0;
  }
  return task;
}

const char *op_idl_vmlmb_msg(int argc, void *argv[])
{
  idl_long_t *size;
  char *csave;

  if (argc != 2) return "expecting 2 arguments";
  size = (idl_long_t *)argv[0];
  if (GET_TYPE(size) != IDL_TYP_DOUBLE)
    return "workspace should be an array of double's";
  if (GET_NELEMENTS(size)*sizeof(double) < OP_MSG_SIZE)
    return "insufficient space in workspace array";
  csave = (char *)argv[1];
  csave[OP_MSG_SIZE] = 0;
  return csave;
}

/*---------------------------------------------------------------------------*/

static char msg_buf[OP_MSG_SIZE];

static void save_error(const char *msg)
{
  int len = OP_MSG_SIZE - 1;
#if 1
  /* skip message prefix */
  if (msg && msg[0] == 'o' && msg[1] == 'p' && msg[2] == '_') {
    int i, c;
    for (i=3 ; i<=len ; ++i) {
      if ((c = msg[i]) == ':') {
	while (++i <= len && msg[i] == ' ')
	  ;
	if (i < len) {
	  msg += i;
	  len -= i;
	}
	break;
      } else if (c == 0) {
	break;
      }
    }
  }
#endif
  memcpy(msg_buf, msg, len);
  msg_buf[len] = 0;
  last_error = msg_buf;
}

static void reset_error()
{
  last_error = 0;
}

static const char *get_error(int reset)
{
  const char *s = last_error;
  if (reset) last_error = 0;
  return s;
}

static void set_error(int *status, const char *message)
{
  if (status) *status = 1;
}

static const char *type_name(idl_long_t type);

static const char *type_name(idl_long_t type)
{
  switch (type) {
  case IDL_TYP_UNDEF: return "UNDEFINED";
  case IDL_TYP_BYTE: return "BYTE";
  case IDL_TYP_FLOAT: return "FLOAT";
  case IDL_TYP_DOUBLE: return "DOUBLE";
  case IDL_TYP_COMPLEX: return "COMPLEX";
  case IDL_TYP_STRING: return "STRING";
  case IDL_TYP_STRUCT: return "STRUCT";
  case IDL_TYP_DCOMPLEX: return "DCOMPLEX";
  case IDL_TYP_PTR: return "POINTER";
  case IDL_TYP_OBJREF: return "OBJECT";
  case IDL_TYP_INT: return "INT";
  case IDL_TYP_UINT: return "UINT";
  case IDL_TYP_LONG: return "LONG";
  case IDL_TYP_ULONG: return "ULONG";
  case IDL_TYP_LONG64: return "LONG64";
  case IDL_TYP_ULONG64: return "ULONG64";
  }
  return "<UNKNOWN_IDL_TYPE>";
}

/*---------------------------------------------------------------------------*/
/* PARSING OF IDL ARGUMENTS (the idea is to use the pair {size(ARG), ARG}
   where ARG is an IDL variable). */

#define FUNC(TYPE, IDL_TYPE_IDENT, IDL_TYPE_NAME)			\
static idl_##TYPE##_t *get_##TYPE##_array(void *size_ptr,		\
					  void *addr,			\
					  idl_long_t *nelements)	\
{									\
  if (! last_error) {							\
    idl_long_t *size = (idl_long_t *)size_ptr;				\
    if (GET_TYPE(size) == IDL_TYPE_IDENT) {				\
      if (nelements) *nelements = GET_NELEMENTS(size);			\
      return (idl_##TYPE##_t *)addr;					\
    }									\
    last_error = "expecting array of " IDL_TYPE_NAME "'s";		\
  }									\
  return 0;								\
}
FUNC(byte,     IDL_TYP_BYTE,     "BYTE")				\
FUNC(int,      IDL_TYP_INT,      "INT")
FUNC(long,     IDL_TYP_LONG,     "LONG")
FUNC(long64,   IDL_TYP_LONG64,   "LONG64")
FUNC(float,    IDL_TYP_FLOAT,    "FLOAT")
FUNC(double,   IDL_TYP_DOUBLE,   "DOUBLE")
FUNC(complex,  IDL_TYP_COMPLEX,  "COMPLEX")
FUNC(dcomplex, IDL_TYP_DCOMPLEX, "DCOMPLEX")

static idl_long_t get_integer(void *size_ptr, void *addr)
{
  if (! last_error) {
    idl_long_t *size = (idl_long_t *)size_ptr;
    if (GET_NDIMS(size) == 0) {
      switch (GET_TYPE(size)) {
      case IDL_TYP_BYTE: return *(idl_byte_t *)addr;
#ifdef idl_int_t
      case IDL_TYP_INT: return *(idl_int_t *)addr;
      case IDL_TYP_UINT: return *(idl_uint_t *)addr;
#endif
#ifdef idl_long_t
      case IDL_TYP_LONG: return *(idl_long_t *)addr;
      case IDL_TYP_ULONG: return *(idl_ulong_t *)addr;
#endif
#ifdef idl_long64_t
      case IDL_TYP_LONG64: return *(idl_long64_t *)addr;
      case IDL_TYP_ULONG64: return *(idl_ulong64_t *)addr;
#endif
      }
    }
    last_error = "expecting scalar integer argument";
  }
  return 0;
}

static idl_double_t get_real(void *size_ptr, void *addr)
{
  if (! last_error) {
    idl_long_t *size = (idl_long_t *)size_ptr;
    if (GET_NDIMS(size) == 0) {
      switch (GET_TYPE(size)) {
      case IDL_TYP_BYTE: return *(idl_byte_t *)addr;
#ifdef idl_int_t
      case IDL_TYP_INT: return *(idl_int_t *)addr;
      case IDL_TYP_UINT: return *(idl_uint_t *)addr;
#endif
#ifdef idl_long_t
      case IDL_TYP_LONG: return *(idl_long_t *)addr;
      case IDL_TYP_ULONG: return *(idl_ulong_t *)addr;
#endif
#ifdef idl_long64_t
      case IDL_TYP_LONG64: return *(idl_long64_t *)addr;
      case IDL_TYP_ULONG64: return *(idl_ulong64_t *)addr;
#endif
      case IDL_TYP_FLOAT: return *(idl_float_t *)addr;
      case IDL_TYP_DOUBLE: return *(idl_double_t *)addr;
      }
    }
    last_error = "expecting scalar real argument";
  }
  return 0;
}

/*---------------------------------------------------------------------------*
 * Local Variables:                                                          *
 * mode: C                                                                   *
 * tab-width: 8                                                              *
 * fill-column: 75                                                           *
 * coding: latin-1                                                           *
 * End:                                                                      *
 *---------------------------------------------------------------------------*/
