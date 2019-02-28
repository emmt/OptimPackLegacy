#
# Makefile --
#
# Makefile for OptimPackLegacy.
#
#------------------------------------------------------------------------------
#
# This file is part of OptimPackLegacy
# <https://github.com/emmt/OptimPackLegacy>.
#
# Copyright (c) 2003-2019, Éric Thiébaut.
#
# OptimPackLegacy is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option) any
# later version.
#
# OptimPackLegacy is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# more details.
#
# You should have received a copy of the GNU General Public License along with
# OptimPackLegacy (file "LICENSE" in the top source directory); if not, write
# to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# 02111-1307 USA
#
#------------------------------------------------------------------------------

srcdir = .

PREFIX =
DESTDIR =
INSTALL = cp -p

CC = gcc -pipe
#CPPFLAGS =  -I. -DOP_INTEGER=long -DOP_LOGICAL=int
CPPFLAGS = -I.
CFLAGS = -O2 -Wall -fPIC

RM = rm -f
AR = ar
ARFLAGS = rv

LIBNAME = liboptimpacklegacy.a
VERSION = `sed < $(srcdir)/../VERSION -e 's/ //g'`
OBJS = opl_algebra.o opl_lnsrch.o opl_utils.o opl_vmlmb.o

all: $(LIBNAME)

install: $(LIBNAME)
	@version=$(VERSION); \
	if test "x$(PREFIX)" = "x"; then \
	  echo "You must define PREFIX macro, e.g.:"; \
	  echo "   > make PREFIX=/usr/local install"; \
	else \
	  test -d "$(PREFIX)/lib" || mkdir -p "$(PREFIX)/lib"; \
	  $(INSTALL) $(LIBNAME) "$(PREFIX)/lib/."; \
	  test -d "$(PREFIX)/include" || mkdir -p "$(PREFIX)/include"; \
	  $(INSTALL) optimpacklegacy.h "$(PREFIX)/include/."; \
	  test -d "$(PREFIX)/doc/OptimPackLegacy-$${version}" || \
	    mkdir -p "$(PREFIX)/doc/OptimPackLegacy-$${version}"; \
	  $(INSTALL) $(srcdir)/../README.md $(srcdir)/../AUTHORS \
	    $(srcdir)/../LICENSE $(srcdir)/optimpacklegacy.h \
	    "$(PREFIX)/doc/OptimPackLegacy-$${version}/."; \
	fi

clean:
	$(RM) *~ $(OBJS) $(LIBNAME)

distclean: clean

$(LIBNAME): $(OBJS)
	$(RM) $(LIBNAME)
	$(AR) $(ARFLAGS) $(LIBNAME) $(OBJS)

%.o: %.c opl_private.h optimpacklegacy.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $(@:.o=.c) -o $@
