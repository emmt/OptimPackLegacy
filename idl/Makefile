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
#	$Id: Makefile,v 1.1 2007/07/05 08:59:03 eric Exp eric $
#	$Log: Makefile,v $
#	Revision 1.1  2007/07/05 08:59:03  eric
#	Initial revision
#
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

OSTYPEVAR = `echo "${OSTYPE}"`
MACHTYPEVAR = `echo "${MACHTYPE}"`

# The following is the default entry point. This section will determine
# what system we are on, set the correct flags and call this same makefile
# again with the correct flags.

all :
	@echo "Building library: $(NAME)$(OSTYPEVAR)"
	@echo "OS type detected: "`uname`
	@echo "Variable OSTYPE : $(OSTYPEVAR)"
	@echo "Variable MACHTYPE : $(MACHTYPEVAR)"
	@case `uname` in \
	"SunOS") if [ \( `/bin/uname -p` != i386 \) \
                      -a \( -x /bin/isainfo \) ]; then \
			rm -f $(OBJS) ;\
			make $(NAME).so \
				"SO_EXT=so" \
				"X_CFLAGS=-fPIC -m64 " \
				"X_LD_FLAGS=-G " ;\
			mv $(NAME).so $(NAME)$(OSTYPEVAR)_64.so ;\
			rm -f $(OBJS) ;\
		fi ;\
		make $(NAME).so \
			"SO_EXT=so" \
			"X_CFLAGS=-fPIC " \
			"X_LD_FLAGS=-G " ;\
			mv $(NAME).so $(NAME)$(OSTYPEVAR).so ;\
		;; \
	"AIX") make $(NAME).a \
		"SO_EXT=a" \
		"CC=cc" \
		"LD=cc" \
		"X_LD_FLAGS=-bM:SRE -bnoentry"\
		"EXPORT_FILES=$(AIX_EXPORT_FILES)" \
		"EXPORT_FLAGS=$(AIX_EXPORT_FLAGS)" ;\
		mv $(NAME).so $(NAME)$(OSTYPEVAR).so ;\
		;; \
	"HP-UX") make $(NAME).sl \
		"SO_EXT=sl" \
		"CC=cc" \
		"LD=ld" \
		"X_CFLAGS=+z -Aa -D_HPUX_SOURCE +e" \
		"X_LD_FLAGS=-b" ;\
		mv $(NAME).so $(NAME)$(OSTYPEVAR).so ;\
		;; \
	"IRIX" ) make $(NAME).so \
		"SO_EXT=so" \
		"CC=cc" \
		"LD=ld" \
		"X_CFLAGS=-n32 -mips3 -KPIC" \
		"X_LD_FLAGS=-shared -n32 -mips3" ;\
		mv $(NAME).so $(NAME)$(OSTYPEVAR).so ;\
		;; \
	"IRIX64" ) make $(NAME).so \
		"SO_EXT=so" \
		"CC=cc" \
		"LD=ld" \
		"X_CFLAGS=-n32 -mips3 -KPIC" \
		"X_LD_FLAGS=-shared -n32 -mips3" ;\
		mv $(NAME).so $(NAME)$(OSTYPEVAR).so ;\
		;; \
	"OSF1" )  make $(NAME).so \
		"SO_EXT=so" \
		"CC=cc" \
		"LD=ld" \
		"X_LD_FLAGS=-S -shared" ;\
		mv $(NAME).so $(NAME)$(OSTYPEVAR).so ;\
		;; \
	"Linux" ) if [ $(MACHTYPEVAR) = x86_64 ]; then \
			rm -f $(OBJS) ;\
			make $(NAME).so \
				"SO_EXT=so" \
				"CC=gcc" \
				"LD=ld" \
				"X_CFLAGS=-fPIC -Wall -m64" \
				"X_LD_FLAGS=-shared" ;\
			mv $(NAME).so $(NAME)$(OSTYPEVAR)_64.so ;\
			rm -f $(OBJS) ;\
		fi ;\
		make $(NAME).so \
		"SO_EXT=so" \
		"CC=gcc" \
		"LD=ld" \
		"X_CFLAGS=-fPIC -Wall" \
		"X_LD_FLAGS=-shared" ;\
		mv $(NAME).so $(NAME)$(OSTYPEVAR).so ;\
		;;\
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
