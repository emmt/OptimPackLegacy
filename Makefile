srcdir = .

PREFIX =
INSTALL = cp -p

CC = gcc -pipe
#CPPFLAGS =  -I. -DOP_INTEGER=long -DOP_LOGICAL=int
CPPFLAGS = -I.
CFLAGS = -O2 -Wall

RM = rm -f
AR = ar
ARFLAGS = rv

LIBNAME = liboptimpack.a
VERSION = `sed < VERSION -e 's/ //g'`
OBJS = op_lnsrch.o op_utils.o op_vmlmb.o
DIRS = yorick idl

DISTRIB_SRC = $(srcdir)
DISTRIB_FILES = AUTHORS LICENSE Makefile NEWS README TODO optimpack.bib

CODE_SRC = $(srcdir)
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

all: $(LIBNAME)
install: $(LIBNAME)
	@version=$(VERSION); \
	if [ "x$(PREFIX)" = "x" ]; then \
	  echo "You must define PREFIX macro, e.g.:"; \
	  echo "   > make PREFIX=/usr/local install"; \
	else \
	  [ -d "$(PREFIX)/lib" ] || mkdir -p "$(PREFIX)/lib"; \
	  $(INSTALL) $(LIBNAME) "$(PREFIX)/lib/."; \
	  [ -d "$(PREFIX)/include" ] || mkdir -p "$(PREFIX)/include"; \
	  $(INSTALL) optimpack.h "$(PREFIX)/include/."; \
	  [ -d "$(PREFIX)/doc/OptimPack-$${version}" ] || \
	    mkdir -p "$(PREFIX)/doc/OptimPack-$${version}"; \
	  $(INSTALL) README AUTHORS LICENSE optimpack.h \
	    "$(PREFIX)/doc/OptimPack-$${version}/."; \
	fi

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
	for file in $(IDL_FILES); do \
	  src=$(IDL_SRC)/$$file; \
	  dst=$$dstdir/$$file; \
	  cp -p "$$src" "$$dst"; \
	  chmod 644 "$$dst"; \
	done; \
	dstdir=$$pkgdir/idl/contrib; \
	mkdir -p "$$dstdir"; \
	for file in $(IDL_CONTRIB_FILES); do \
	  src=$(IDL_CONTRIB_SRC)/$$file; \
	  dst=$$dstdir/$$file; \
	  cp -p "$$src" "$$dst"; \
	  chmod 644 "$$dst"; \
	done; \
	dstdir=$$pkgdir/yorick; \
	mkdir -p "$$dstdir"; \
	for file in $(YORICK_FILES); do \
	  src=$(YORICK_SRC)/$$file; \
	  dst=$$dstdir/$$file; \
	  cp -p "$$src" "$$dst"; \
	  chmod 644 "$$dst"; \
	done; \
	dstdir=$$pkgdir; \
	mkdir -p "$$dstdir"; \
	for file in $(DISTRIB_FILES); do \
	  src=$(DISTRIB_SRC)/$$file; \
	  dst=$$dstdir/$$file; \
	  cp -p "$$src" "$$dst"; \
	  chmod 644 "$$dst"; \
	done; \
	dst=$$dstdir/VERSION; \
	echo "$$version" >"$$dst"; \
	chmod 644 "$$dst"; \
	dstdir=$$pkgdir; \
	mkdir -p "$$dstdir"; \
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
	$(RM) *~ $(OBJS) $(LIBNAME)
	for dir in $(DIRS); do \
	  if [ -f "$$dir/Makefile" ]; then \
	    (cd $$dir; make clean); \
	  fi; \
	done

distclean: clean
	$(RM) $(LIBNAME)

$(LIBNAME): $(OBJS)
	$(RM) $(LIBNAME)
	$(AR) $(ARFLAGS) $(LIBNAME) $(OBJS)

op_lnsrch.o: op_lnsrch.c optimpack.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $(@:.o=.c) -o $@
op_vmlmb.o: op_vmlmb.c optimpack.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $(@:.o=.c) -o $@
op_cgmnb.o: op_cgmnb.c optimpack.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $(@:.o=.c) -o $@
op_utils.o: op_utils.c optimpack.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $(@:.o=.c) -o $@
