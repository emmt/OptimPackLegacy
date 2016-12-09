/*
 * op_limits.h --
 *
 *	Definitions to guess sizes of integers on this machine.
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

#ifndef _OP_LIMITS_H
#define _OP_LIMITS_H 1

#include <limits.h>

#define OP_INT16_MIN               (-32768)
#define OP_INT16_MAX                 32767
#define OP_UINT16_MAX                65535

#define OP_INT32_MIN          (-2147483648L)
#define OP_INT32_MAX            2147483647L
#define OP_UINT32_MAX           4294967295U

#define OP_INT64_MIN (-9223372036854775808L)
#define OP_INT64_MAX   9223372036854775807L
#define OP_UINT64_MAX 18446744073709551615UL

#if ! defined (OP_SHRT_BITS) && defined(USHRT_MAX)
# if (USHRT_MAX == OP_UINT16_MAX)
#  define OP_SHRT_BITS    16
# else
#  if (USHRT_MAX == OP_UINT32_MAX)
#   define OP_SHRT_BITS   32
#  else
#   if (USHRT_MAX == OP_UINT64_MAX)
#    define OP_SHRT_BITS  64
#   endif
#  endif
# endif
#endif
#if ! defined (OP_SHRT_BITS) && defined(SHRT_MIN) && defined(SHRT_MAX)
# if (SHRT_MIN == OP_INT16_MIN) && SHRT_MAX == OP_INT16_MAX)
#  define OP_SHRT_BITS    16
# else
#  if (SHRT_MIN == OP_INT32_MIN) && SHRT_MAX == OP_INT32_MAX)
#   define OP_SHRT_BITS   32
#  else
#   if (SHRT_MIN == OP_INT64_MIN) && SHRT_MAX == OP_INT64_MAX)
#    define OP_SHRT_BITS  64
#   endif
#  endif
# endif
#endif

#if ! defined (OP_INT_BITS) && defined(UINT_MAX)
# if (UINT_MAX == OP_UINT16_MAX)
#  define OP_INT_BITS    16
# else
#  if (UINT_MAX == OP_UINT32_MAX)
#   define OP_INT_BITS   32
#  else
#   if (UINT_MAX == OP_UINT64_MAX)
#    define OP_INT_BITS  64
#   endif
#  endif
# endif
#endif
#if ! defined (OP_INT_BITS) && defined(INT_MIN) && defined(INT_MAX)
# if (INT_MIN == OP_INT16_MIN) && (INT_MAX == OP_INT16_MAX)
#  define OP_INT_BITS    16
# else
#  if (INT_MIN == OP_INT32_MIN) && (INT_MAX == OP_INT32_MAX)
#   define OP_INT_BITS   32
#  else
#   if (INT_MIN == OP_INT64_MIN) && (INT_MAX == OP_INT64_MAX)
#    define OP_INT_BITS  64
#   endif
#  endif
# endif
#endif

#if ! defined (OP_LONG_BITS) && defined(ULONG_MAX)
# if (ULONG_MAX == OP_UINT16_MAX)
#  define OP_LONG_BITS    16
# else
#  if (ULONG_MAX == OP_UINT32_MAX)
#   define OP_LONG_BITS   32
#  else
#   if (ULONG_MAX == OP_UINT64_MAX)
#    define OP_LONG_BITS  64
#   endif
#  endif
# endif
#endif
#if ! defined (OP_LONG_BITS) && defined(LONG_MIN) && defined(LONG_MAX)
# if (LONG_MIN == OP_INT16_MIN) && (LONG_MAX == OP_INT16_MAX)
#  define OP_LONG_BITS    16
# else
#  if (LONG_MIN == OP_INT32_MIN) && (LONG_MAX == OP_INT32_MAX)
#   define OP_LONG_BITS   32
#  else
#   if (LONG_MIN == OP_INT64_MIN) && (LONG_MAX == OP_INT64_MAX)
#    define OP_LONG_BITS  64
#   endif
#  endif
# endif
#endif

#if ! defined (OP_LLONG_BITS) && defined(ULLONG_MAX)
# if (ULLONG_MAX == OP_UINT16_MAX)
#  define OP_LLONG_BITS    16
# else
#  if (ULLONG_MAX == OP_UINT32_MAX)
#   define OP_LLONG_BITS   32
#  else
#   if (ULLONG_MAX == OP_UINT64_MAX)
#    define OP_LLONG_BITS  64
#   endif
#  endif
# endif
#endif
#if ! defined (OP_LLONG_BITS) && defined(LLONG_MIN) && defined(LLONG_MAX)
# if (LLONG_MIN == OP_INT16_MIN) && (LLONG_MAX == OP_INT16_MAX)
#  define OP_LLONG_BITS    16
# else
#  if (LLONG_MIN == OP_INT32_MIN) && (LLONG_MAX == OP_INT32_MAX)
#   define OP_LLONG_BITS   32
#  else
#   if (LLONG_MIN == OP_INT64_MIN) && (LLONG_MAX == OP_INT64_MAX)
#    define OP_LLONG_BITS  64
#   endif
#  endif
# endif
#endif

/*---------------------------------------------------------------------------*/
#endif /* _OP_LIMITS_H */
