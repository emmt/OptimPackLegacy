#
# Makefile --
#
# Makefile for OptimPack.
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
DISTRIB_FILES = AUTHORS LICENSE Makefile NEWS.md README.md TODO.md optimpack.bib

CODE_SRC = $(srcdir)/src
CODE_FILES = op_limits.h op_lnsrch.c optimpack.h op_utils.c op_vmlmb.c

IDL_SRC = $(srcdir)/idl
IDL_FILES = \
    Makefile \
    op_ensure_double.pro \
    op_ensure_long.pro \
    op_idl.c \
    op_idl.h \
    op_init.pro \
    op_last_error.pro \
    op_test.pro \
    op_vmlmb_msg.pro \
    op_vmlmb.pro \
    op_vmlmb_setup.pro \
    op_wrapper.c \
    README

IDL_CONTRIB_SRC = $(srcdir)/idl/contrib
IDL_CONTRIB_FILES = fmin_op.pro  Makefile.OptimPack

YORICK_SRC = $(srcdir)/yorick
YORICK_FILES = Makefile \
    lbfgsb.f lbfgsb.i lbfgsb_wrapper.c \
    lbfgs.f lbfgs.i lbfgs_wrapper.c \
    OptimPack1.i OptimPack1-test.i optimpack.c

all:
	@echo "No default target"

distrib:
	@version=$(VERSION); \
	if test "x$$version" = "x"; then \
	  echo >&2 "bad VERSION"; \
	  return 1; \
	fi; \
	pkgdir=optimpack-$$version; \
	archive=$$pkgdir.tar.bz2; \
	if test -e "$$pkgdir"; then \
	  echo >&2 "error: $$pkgdir already exists"; \
	  return 1; \
	fi; \
	if test -e "$$archive"; then \
	  echo >&2 "error: $$archive already exists"; \
	  return 1; \
	fi; \
	dstdir=$$pkgdir/idl; \
	mkdir -p "$$dstdir"; \
	chmod 755 "$$dstdir"; \
	for file in $(IDL_FILES); do \
	  src=$(IDL_SRC)/$$file; \
	  dst=$$dstdir/$$file; \
	  cp -p "$$src" "$$dst"; \
	  chmod 644 "$$dst"; \
	done; \
	dstdir=$$pkgdir/idl/contrib; \
	mkdir -p "$$dstdir"; \
	chmod 755 "$$dstdir"; \
	for file in $(IDL_CONTRIB_FILES); do \
	  src=$(IDL_CONTRIB_SRC)/$$file; \
	  dst=$$dstdir/$$file; \
	  cp -p "$$src" "$$dst"; \
	  chmod 644 "$$dst"; \
	done; \
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

