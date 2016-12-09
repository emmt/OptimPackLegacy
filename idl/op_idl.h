/*
 * op_idl.h --
 *
 *	Definitions for OptimPack-IDL interface.
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

#ifndef _OP_IDL_H
#define _OP_IDL_H 1


/* Definition of types matching the ones in IDL:
    C-type         IDL type       IDL code
    idl_int_t        INT        IDL_TYP_INT
    idl_long_t       LONG       IDL_TYP_LONG
    idl_long64_t     LONG64     IDL_TYP_LONG64

*/
#include "stdio.h" /* needed by IDL's export.h */
#include "export.h"
#define idl_byte_t     unsigned char
#define idl_int_t      short
#define idl_uint_t     unsigned idl_int_t
#define idl_long_t     IDL_LONG
#define idl_ulong_t    IDL_ULONG
#define idl_long64_t   IDL_LONG64
#define idl_ulong64_t  IDL_ULONG64
#define idl_float_t    float
#define idl_double_t   double
#define idl_complex_t  IDL_COMPLEX
#define idl_dcomplex_t IDL_DCOMPLEX

/*---------------------------------------------------------------------------*/

/*
 * Include "optimpack.h" _after_ having set the definitions of the
 * various types used by OptimPack routines.
 */

#define OP_INTEGER idl_long_t
#define OP_LOGICAL idl_byte_t


/*---------------------------------------------------------------------------*/
extern const char *op_idl_last_error(int argc, void *argv[]);
extern idl_long_t  op_idl_vmlmb_first(int argc, void *argv[]);
extern idl_long_t  op_idl_vmlmb_next(int argc, void *argv[]);
extern const char *op_idl_vmlmb_msg(int argc, void *argv[]);



/*---------------------------------------------------------------------------*/
#endif /* _OP_IDL_H */
