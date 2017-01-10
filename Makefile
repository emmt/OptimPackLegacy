#
# Makefile --
#
# Makefile for OptimPackLegacy.
#
#------------------------------------------------------------------------------
#
# Copyright (c) 2003, 2016 Éric Thiébaut.
#
# This file is part of OptimPack <https://github.com/emmt/OptimPackLegacy>.
#
# OptimPack is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# OptimPack is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# more details.
#
# You should have received a copy of the GNU General Public License along
# with OptimPack (file "LICENSE" in the top source directory); if not,
# write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
#
#------------------------------------------------------------------------------

srcdir = .

VERSION = `sed < VERSION -e 's/ //g'`
SUBDIRS = yorick idl src

DISTRIB_SRC = $(srcdir)
DISTRIB_FILES = AUTHORS LICENSE Makefile NEWS.md README.md TODO.md \
    optimpack.bib

CODE_SRC = $(srcdir)/src
CODE_FILES = opl_algebra.c opl_lnsrch.c opl_vmlmb.c opl_utils.c \
     opl_limits.h  opl_private.h optimpacklegacy.h

YORICK_SRC = $(srcdir)/yorick
YORICK_FILES = Makefile opl_yorick.c optimpacklegacy.i \
    optimpacklegacy-start.i \
    optimpacklegacy-tests.i optimpacklegacy-tests.out

all:
	@echo "No default target"

distrib:
	@version=$(VERSION); \
	if test "x$$version" = "x"; then \
	  echo >&2 "bad VERSION"; \
	  return 1; \
	fi; \
	pkgdir=OptimPackLegacy-$$version; \
	archive=$$pkgdir.tar.bz2; \
	if test -e "$$pkgdir"; then \
	  echo >&2 "error: $$pkgdir already exists"; \
	  return 1; \
	fi; \
	if test -e "$$archive"; then \
	  echo >&2 "error: $$archive already exists"; \
	  return 1; \
	fi; \
	dstdir=$$pkgdir/yorick; \
	mkdir -p "$$dstdir"; \
	chmod 755 "$$dstdir"; \
	for file in $(YORICK_FILES); do \
	  src=$(YORICK_SRC)/$$file; \
	  dst=$$dstdir/$$file; \
	  cp -p "$$src" "$$dst"; \
	  chmod 644 "$$dst"; \
	done; \
	dstdir=$$pkgdir; \
	mkdir -p "$$dstdir"; \
	chmod 755 "$$dstdir"; \
	for file in $(DISTRIB_FILES); do \
	  src=$(DISTRIB_SRC)/$$file; \
	  dst=$$dstdir/$$file; \
	  cp -p "$$src" "$$dst"; \
	  chmod 644 "$$dst"; \
	done; \
	dst=$$dstdir/VERSION; \
	echo "$$version" >"$$dst"; \
	chmod 644 "$$dst"; \
	dstdir=$$pkgdir/src; \
	mkdir -p "$$dstdir"; \
	chmod 755 "$$dstdir"; \
	for file in $(CODE_FILES); do \
	  src=$(CODE_SRC)/$$file; \
	  dst=$$dstdir/$$file; \
	  cp -p "$$src" "$$dst"; \
	  chmod 644 "$$dst"; \
	done; \
	tar cf - "$$pkgdir" | bzip2 -9 > "$$archive"; \
	rm -rf "$$pkgdir"; \
	echo "archive $$archive created"; \
	return 0

clean:
	$(RM) *~
	for dir in $(SUBDIRS); do \
	  if [ -f "$$dir/Makefile" ]; then \
	    (cd $$dir; make clean); \
	  fi; \
	done

distclean: clean
	for dir in $(SUBDIRS); do \
	  if [ -f "$$dir/Makefile" ]; then \
	    (cd $$dir; make distclean); \
	  fi; \
	done

