/*
 * opl_utils.c --
 *
 * Utilities routines for OptimPackLegacy library.
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

#include <string.h>
#include <stdio.h>
#include "optimpacklegacy.h"

/*---------------------------------------------------------------------------*/

int opl_error(char *buf, const char *errmsg)
{
  if (buf) {
    if (! errmsg) errmsg = "unknown error";
    strncpy(buf, errmsg, OPL_MSG_LEN);
    buf[OPL_MSG_LEN] = 0;
  }
  return OPL_ERROR;
}

void opl_mcopy(const char *msg, char *buf)
{
  if (buf) {
    if (msg) {
      strncpy(buf, msg, OPL_MSG_LEN);
      buf[OPL_MSG_LEN] = 0;
    } else {
      buf[0] = 0;
    }
  }
}

/*---------------------------------------------------------------------------*/
