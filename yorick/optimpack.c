/*
 * optimpack.c --
 *
 *	Implementation of Yorick wrapper for OptimPack.
 *
 *-----------------------------------------------------------------------------
 *
 *	Copyright (C) 2003-2007 Eric Thiébaut.
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
 * History:
 *	$Id$
 *	$Log$
 *-----------------------------------------------------------------------------
 */

void Y_op_vmlmb_first(int argc)
{
  Dimension *dims;
  Symbol *argv = sp - argc + 1;
  long n, m;
  double fmin, fatol, frtol, sftol, sgtol, sxtol;
  char *csave;
  long *isave;
  double *dsave;
  long ncsave, nisave, ndsave;
  long task;
  if (argc != 11) YError("op_vmlmb_first takes exactly 11 arguments");
  n = YGetInteger(&argv[0]);
  m = YGetInteger(&argv[1]);
  fmin = YGetReal(&argv[2]);
  fatol = YGetReal(&argv[3]);
  frtol = YGetReal(&argv[4]);
  sftol = YGetReal(&argv[5]);
  sgtol = YGetReal(&argv[6]);
  sxtol = YGetReal(&argv[7]);
  csave = get_char_array(&argv[8], 0, &ncsave);
  isave = get_long_array(&argv[9], 0, &nisave);
  dsave = get_double_array(&argv[10], 0, &ndsave);
  if (ncsave < OP_MSG_SIZE) YError("too few elements for CSAVE");
  if (nisave < ???) YError("too few elements for ISAVE");
  if (ndsave < ???) YError("too few elements for DSAVE");
  task = op_vmlmb_first(n, m, fmin, fatol, frtol, sftol, sgtol, sxtol,
			csave, isave, dsave);
  if (task != OP_TASK_FG) YError(csave);
  PushLongValue(task);
}

static void parse_ws(Symbol *s,
		     char   **csave_ptr, long ncsave,
		     long   **isave_ptr, long nisave,
		     double **dsave_ptr, long ndsave)
{
  void *ptr;
  Operand op;
  if (!s->ops) YError("unexpected keyword argument");
  s->ops->FormOperand(s, &op);
  if (op.ops->typeID!=T_POINTER || op.type.number != 3)
    YError("expecting array of 3 pointers");
  ptr = *(void **)op.value;

}

#define IS_REF(S)       ((S)->ops == &referenceSym)
#define DEREF_SYMBOL(S) (IS_REF(S) ? &globTab[(S)->index] : (S))

static void unexpected_keyword_argument(const char *name)
{
  if (name && strlen(name) < 30) {
    char msg[80];
    strcpy(msg, "unexpected keyword argument for ");
    strcat(msg, name);
    YError(msg);
  } else {
    YError("unexpected keyword argument");
  }
}

static char *get_char_array(const char *name, Symbol *s,
			    int nil_ok, long *number);
static char *get_char_array(const char *name, Symbol *s,
			    int nil_ok, long *number)
{
  Operand op;
  if (! s->ops) YError("unexpected keyword argument");
  s->ops->FormOperand(s, &op);
  if (op.ops->typeID == T_CHAR) {
    if (nil_ok && ) {
      if (number) *number = 0;
      return NULL;
    }
    YError("expecting character array argument");
  }
  if (number) *number = op.type.number;
  return ;
  if (op.ops->typeID!=T_POINTER || op.type.dims)
    YError("expecting scalar pointer argument");
  return *(void **)op.value;
}

/*---------------------------------------------------------------------------*
 * Local Variables:                                                          *
 * mode: C                                                                   *
 * tab-width: 8                                                              *
 * fill-column: 75                                                           *
 * coding: latin-1                                                           *
 * End:                                                                      *
 *---------------------------------------------------------------------------*/
