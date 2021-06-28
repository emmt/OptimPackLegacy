/*
 * opl_utils.c --
 *
 * Utilities routines for OptimPackLegacy library.
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

#include <stdarg.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include "opl_private.h"

/*---------------------------------------------------------------------------*/
/* CONTEXT MANAGEMENT AND ERROR REPORTING */

const char* _opl_success_message = "Successful operation";
static const char* _opl_insufficient_memory_message = "Insufficient memory";
static const char* _opl_illegal_address_message = "Illegal address";
static const char* _opl_invalid_argument_message = "Invalid argument";
static const char* _opl_out_of_bounds_message = "Out of bounds size or index";
static const char* _opl_corrupted_message = "Corrupted data";
static const char* _opl_overflow_message = "Numerical overflow";
#if 0
static const char* _opl_syntax_error_message = "Syntax error";
static const char* _opl_open_error_message = "Cannot open file";
static const char* _opl_close_error_message = "Error while closing file";
static const char* _opl_read_error_message = "Read error";
static const char* _opl_write_error_message = "Write error";
static const char* _opl_stat_error_message = "Cannot stat file";
static const char* _opl_file_exists_message = "File already exists";
static const char* _opl_zip_error_message = "Error with PkZip archive";
#endif

extern void opl_initialize_context(
    opl_context_t* ctx)
{
    ctx->message = _opl_success_message;
    ctx->status = OPL_SUCCESS;
    ctx->syserr = 0;
}

#define DEFINE_FUNCTION(name, NAME)                     \
    opl_status_t opl_##name(                            \
        opl_context_t* ctx)                             \
    {                                                   \
        return _OPL_SET_CONTEXT(                        \
            ctx, OPL_##NAME, _opl_##name##_message);    \
    }
DEFINE_FUNCTION(success,             SUCCESS)
DEFINE_FUNCTION(insufficient_memory, INSUFFICIENT_MEMORY)
DEFINE_FUNCTION(illegal_address,     ILLEGAL_ADDRESS)
DEFINE_FUNCTION(invalid_argument,    INVALID_ARGUMENT)
DEFINE_FUNCTION(out_of_bounds,       OUT_OF_BOUNDS)
DEFINE_FUNCTION(corrupted,           CORRUPTED)
DEFINE_FUNCTION(overflow,            OVERFLOW)
#if 0
DEFINE_FUNCTION(syntax_error,        SYNTAX_ERROR)
DEFINE_FUNCTION(open_error,          OPEN_ERROR)
DEFINE_FUNCTION(close_error,         CLOSE_ERROR)
DEFINE_FUNCTION(read_error,          READ_ERROR)
DEFINE_FUNCTION(write_error,         WRITE_ERROR)
DEFINE_FUNCTION(stat_error,          STAT_ERROR)
DEFINE_FUNCTION(file_exists,         FILE_EXISTS)
DEFINE_FUNCTION(zip_error,           ZIP_ERROR)
#endif
#undef DEFINE_FUNCTION

    opl_status_t opl_system_error(
        opl_context_t* ctx)
{
    int syserr = errno;
    ctx->syserr = syserr;
    ctx->message = strerror(syserr);
    return (ctx->status = OPL_SYSTEM_ERROR);
}

opl_status_t opl_get_status(
    opl_context_t* ctx)
{
    return ctx->status;
}

int opl_get_errno(
    opl_context_t* ctx)
{
    return (ctx->status == OPL_SYSTEM_ERROR ? ctx->syserr : 0);
}

const char* opl_get_message(
    opl_context_t* ctx)
{
    return ctx->message;
}

#define TRUNCATE(str, siz)                      \
    (str)[(siz) - 6] = '[';                     \
    (str)[(siz) - 5] = '.';                     \
    (str)[(siz) - 4] = '.';                     \
    (str)[(siz) - 3] = '.';                     \
    (str)[(siz) - 2] = ']';                     \
    (str)[(siz) - 1] = '\0'

opl_status_t opl_set_context(
    opl_context_t* ctx,
    opl_status_t status,
    const char* message,
    opl_storage_type_t type)
{
    if (type == OPL_PERMANENT) {
        ctx->message = message;
    } else {
        size_t size = sizeof(ctx->buffer);
        size_t nchars = (message == NULL ? 0 : strlen(message));
        if (nchars <= 0) {
            ctx->buffer[0] = '\0';
        } else if (nchars < size) {
            memcpy(ctx->buffer, message, nchars + 1);
        } else {
            memcpy(ctx->buffer, message, size);
            TRUNCATE(ctx->buffer, size);
        }
        ctx->message = ctx->buffer;
    }
    ctx->syserr = (status == OPL_SYSTEM_ERROR ? errno : 0);
    return (ctx->status = status);
}

opl_status_t opl_format_context(
    opl_context_t* ctx,
    opl_status_t status,
    const char* format,
    ...)
{
    const size_t size = sizeof(ctx->buffer);
    va_list ap;
    size_t nchars;

    va_start(ap, format);
    nchars = vsnprintf(ctx->buffer, size, format, ap);
    va_end(ap);
    if (nchars >= size) {
        TRUNCATE(ctx->buffer, size);
    }
    ctx->message = ctx->buffer;
    ctx->syserr = (status == OPL_SYSTEM_ERROR ? errno : 0);
    return (ctx->status = status);
}

const char* opl_get_default_message(
    opl_status_t status)
{
    switch(status) {
    case OPL_SUCCESS:             return _opl_success_message;
    case OPL_INSUFFICIENT_MEMORY: return _opl_insufficient_memory_message;
    case OPL_ILLEGAL_ADDRESS:     return _opl_illegal_address_message;
    case OPL_INVALID_ARGUMENT:    return _opl_invalid_argument_message;
    case OPL_OUT_OF_BOUNDS:       return _opl_out_of_bounds_message;
    case OPL_CORRUPTED:           return _opl_corrupted_message;
    case OPL_OVERFLOW:            return _opl_overflow_message;
#if 0
    case OPL_SYNTAX_ERROR:        return _opl_syntax_error_message;
    case OPL_OPEN_ERROR:          return _opl_open_error_message;
    case OPL_CLOSE_ERROR:         return _opl_close_error_message;
    case OPL_READ_ERROR:          return _opl_read_error_message;
    case OPL_WRITE_ERROR:         return _opl_write_error_message;
    case OPL_STAT_ERROR:          return _opl_stat_error_message;
    case OPL_FILE_EXISTS:         return _opl_file_exists_message;
    case OPL_ZIP_ERROR:           return _opl_zip_error_message;
#endif
    default:                      return "Unknown status";
    }
}

/*---------------------------------------------------------------------------*/
