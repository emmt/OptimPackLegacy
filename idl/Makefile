#
# Makefile -
#
#	Makefile for OptimPack / IDL plugin.
#
#------------------------------------------------------------------------------
#
#	Copyright (c) 1988-2001, Research Systems Inc.  All rights reserved.
#	This software includes information which is proprietary to and a
#	trade secret of Research Systems, Inc.  It is not to be disclosed
#	to anyone outside of this organization. Reproduction by any means
#	whatsoever is  prohibited without express written permission.
#
#	makefile for CALL_EXTERNAL examples, UNIX platforms
#
#	This is a recursive makefile. The default target "all" detects which
#	operating system you are on, sets system specific macros, and then
#	calls make again to build the target "libs".  Specifying a target
#	other than "all","clean", or "tidy" will NOT set all required macros
#	and will probably produce incorrect results.
#
#------------------------------------------------------------------------------
#
#	$Id: makefile,v 1.9 2001/01/15 22:27:36 scottm Exp $
#	$Log$
#------------------------------------------------------------------------------

SHELL=/bin/sh
# 	The following macro might need to be edited if the location
#	of IDL is not the standard location on your system.
#IDL_DIR = /usr/local/rsi/idl_5.5


# The following macros establish general defaults. They are overridden
# in the "all" rule as needed to satisfy a given platform's requirements.

CC		= gcc
CFLAGS		= -O2
OP_CFLAGS	= -I.. -I. -I$(IDL_DIR)/external $(CFLAGS) $(X_CFLAGS)
LD 		= ld
SHELL		= /bin/sh
#platform specific C compiler flags (redefined below)
X_CFLAGS	=
#platform specific linker flags
X_LD_FLAGS	=

SO_EXT		=so

NAME		=OptimPack_IDL



# The following is the default entry point. This section will determine
# what system we are on, set the correct flags and call this same makefile
# again with the correct flags.

all :
	@echo "OS type detected: "`uname`
	@case `uname` in \
	"SunOS") if [ \( `/bin/uname -p` != i386 \) \
                      -a \( -x /bin/isainfo \) ]; then \
			rm -f $(OBJS) ;\
			make $(NAME).so \
				"SO_EXT=so" \
				"CC=cc -xtarget=ultra -xarch=v9" \
				"LD=cc -xtarget=ultra -xarch=v9" \
				"X_CFLAGS=-K pic " \
				"X_LD_FLAGS=-G " ;\
			mv $(NAME).so $(NAME)_64.so ;\
			rm -f $(OBJS) ;\
		fi ;\
		make $(NAME).so \
			"SO_EXT=so" \
			"CC=cc" \
			"LD=cc" \
			"X_CFLAGS=-K pic " \
			"X_LD_FLAGS=-G " ;\
		;; \
	"AIX") make $(NAME).a \
		"SO_EXT=a" \
		"CC=cc" \
		"LD=cc" \
		"X_LD_FLAGS=-bM:SRE -bnoentry"\
		"EXPORT_FILES=$(AIX_EXPORT_FILES)" \
		"EXPORT_FLAGS=$(AIX_EXPORT_FLAGS)";; \
	"HP-UX") make $(NAME).sl \
		"SO_EXT=sl" \
		"CC=cc" \
		"LD=ld" \
		"X_CFLAGS=+z -Aa -D_HPUX_SOURCE +e" \
		"X_LD_FLAGS=-b";; \
	"IRIX" ) make $(NAME).so \
		"SO_EXT=so" \
		"CC=cc" \
		"LD=ld" \
		"X_CFLAGS=-n32 -mips3 -KPIC" \
		"X_LD_FLAGS=-shared -n32 -mips3";; \
	"IRIX64" ) make $(NAME).so \
		"SO_EXT=so" \
		"CC=cc" \
		"LD=ld" \
		"X_CFLAGS=-n32 -mips3 -KPIC" \
		"X_LD_FLAGS=-shared -n32 -mips3";; \
	"OSF1" )  make $(NAME).so \
		"SO_EXT=so" \
		"CC=cc" \
		"LD=ld" \
		"X_LD_FLAGS=-S -shared";; \
	"Linux" ) make $(NAME).so \
		"SO_EXT=so" \
		"CC=gcc" \
		"LD=ld" \
		"X_CFLAGS=-fPIC -Wall" \
		"X_LD_FLAGS=-shared" ;;\
	*) echo "This system is not supported" ;; \
	esac


OBJS = op_idl.o op_utils.o op_vmlmb.o op_lnsrch.o

AIX_EXPORT_FLAGS= \
  -bE:op_idl_last_error.export \
  -bE:op_idl_vmlmb_first.export \
  -bE:op_idl_vmlmb_next.export \
  -bE:op_idl_vmlmb_msg.export

AIX_EXPORT_FILES= \
  op_idl_last_error.export \
  op_idl_vmlmb_first.export \
  op_idl_vmlmb_next.export \
  op_idl_vmlmb_msg.export

op_idl_last_error.export:
	echo $(@:.export=) > $(@)
op_idl_vmlmb_first.export:
	echo $(@:.export=) > $(@)
op_idl_vmlmb_msg.export:
	echo $(@:.export=) > $(@)
op_idl_vmlmb_next.export:
	echo $(@:.export=) > $(@)

# Rules for object files:
.c.o : op_idl.h
	$(CC) $(OP_CFLAGS) -c $*.c

# Rules for OptimPack object files:
op_lnsrch.o: ../op_lnsrch.c ../optimpack.h op_wrapper.c op_idl.h
	$(CC) $(OP_CFLAGS) -c op_wrapper.c -DSOURCE=\"../$(@:.o=.c)\" -o $@
op_vmlmb.o: ../op_vmlmb.c ../optimpack.h op_wrapper.c op_idl.h
	$(CC) $(OP_CFLAGS) -c op_wrapper.c -DSOURCE=\"../$(@:.o=.c)\" -o $@
op_cgmnb.o: ../op_cgmnb.c ../optimpack.h op_wrapper.c op_idl.h
	$(CC) $(OP_CFLAGS) -c op_wrapper.c -DSOURCE=\"../$(@:.o=.c)\" -o $@
op_utils.o: ../op_utils.c ../optimpack.h op_wrapper.c op_idl.h
	$(CC) $(OP_CFLAGS) -c op_wrapper.c -DSOURCE=\"../$(@:.o=.c)\" -o $@

$(NAME).$(SO_EXT): $(OBJS) $(EXPORT_FILES)
	$(LD) $(X_LD_FLAGS) -o $(NAME).$(SO_EXT) $(OBJS) $(EXPORT_FLAGS)

tidy:
	rm -f *.o *~

clean: tidy
	rm -f *.export *.so *.sl *.a
