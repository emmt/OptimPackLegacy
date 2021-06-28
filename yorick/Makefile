#
# Makefile --
#
# Makefile for YOPL (Yorick + OptimPackLegacy).
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

# ---------------------------------------------------------------- Yorick setup
# these values filled in by `yorick -batch make.i`
Y_MAKEDIR=
Y_EXE=
Y_EXE_PKGS=
Y_EXE_HOME=
Y_EXE_SITE=
Y_HOME_PKG=

# ---------------------------------------------------------- optimization flags

# options for make command line, e.g. `make COPT=-g TGT=exe`
COPT = $(COPT_DEFAULT)
TGT = $(DEFAULT_TGT)

# ----------------------------------------------------- macros for this package

PKG_NAME = optimpacklegacy
PKG_I = $(srcdir)/optimpacklegacy.i

OBJS = opl_algebra.o opl_lnsrch.o opl_utils.o opl_vmlmb.o opl_yorick.o

# change to give the executable a name other than yorick
PKG_EXENAME = yorick

# `PKG_DEPLIBS=-Lsomedir -lsomelib` for dependencies of this package
PKG_DEPLIBS =
# set compiler (or rarely loader) flags specific to this package
PKG_CFLAGS =
PKG_LDFLAGS =

# list of additional package names you want in PKG_EXENAME
# (typically Y_EXE_PKGS should be first here)
EXTRA_PKGS = $(Y_EXE_PKGS)

# list of additional files for clean
PKG_CLEAN =

# autoload file for this package, if any
PKG_I_START = $(srcdir)/optimpacklegacy-start.i

# non-pkg.i include files for this package, if any
PKG_I_EXTRA =

# ------------------------------------- standard targets and rules (in Makepkg)

# set macros Makepkg uses in target and dependency names
# DLL_TARGETS, LIB_TARGETS, EXE_TARGETS
# are any additional targets (defined below) prerequisite to
# the plugin library, archive library, and executable, respectively
PKG_I_DEPS = $(PKG_I)
Y_DISTMAKE = distmake

ifeq (,$(strip $(Y_MAKEDIR)))
$(info *** WARNING: Y_MAKEDIR not defined, you may run 'yorick -batch make.i' first)
else
include $(Y_MAKEDIR)/Make.cfg
include $(Y_MAKEDIR)/Makepkg
include $(Y_MAKEDIR)/Make$(TGT)
endif

# override macros Makepkg sets for rules and other macros
# Y_HOME and Y_SITE in Make.cfg may not be correct (e.g.- relocatable)
Y_HOME = $(Y_EXE_HOME)
Y_SITE = $(Y_EXE_SITE)

# reduce chance of yorick-1.5 corrupting this Makefile
MAKE_TEMPLATE = protect-against-1.5

# ------------------------------------------ targets and rules for this package

# Directory where are the OptimPackLegacy sources:
OPT_SRC = $(srcdir)/../src

# simple example:
#myfunc.o: myapi.h
# more complex example (also consider using PKG_CFLAGS above):
#myfunc.o: myapi.h myfunc.c
#	$(CC) $(CPPFLAGS) $(CFLAGS) -DMY_SWITCH -o $@ -c myfunc.c

opl_yorick.o: $(srcdir)/opl_yorick.c $(OPT_SRC)/optimpacklegacy.h
	$(CC) $(CPPFLAGS) -I$(OPT_SRC) $(CFLAGS) -o $@ -c $<

opl_%.o: $(OPT_SRC)/opl_%.c $(OPT_SRC)/opl_private.h $(OPT_SRC)/optimpacklegacy.h
	$(CC) $(CPPFLAGS) -I$(OPT_SRC) $(CFLAGS) -o $@ -c $<

# ------------------------------------------------------------- end of Makefile
