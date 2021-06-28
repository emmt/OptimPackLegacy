# OptimPackLegacy (version 1)

This is **OptimPackLegacy**, a C library for optimization of large scale
problems possibly with bound constraints.  This version implements:

- Moré and Thuente method for inexact line search;

- VMLMB algorithm by Éric Thiébaut which is a limited memory BFGS (variable
  metric) method possibly with bound constraints and/or preconditioning.

In order to make embedding OptimPackLegacy into another language as easy as
possible, the routines use reverse communication: all local variables needed by
the optimization routines get saved into workspace structures (allocated by the
library or provided by the caller) and the optimization routines never
explicitely call the penalty function to optimize.

Most of the documention is in the header file
[optimpacklegacy.h](src/optimpacklegacy.h).

Another version of OptimPack is under development and available
[here](https://github.com/emmt/OptimPack).  This new version implements a
completely different memory management system to be more flexible and allow for
various kinds of storage (standard memory, distributed memory, GPU, *etc.*).

Directory [yorick](yorick) provides a Yorick interface to OptimPackLegacy.

Directory [python](python) provides a Python interface to OptimPackLegacy.


## References

* Jorge J. Moré and David J. Thuente, "*Line search algorithms with guaranteed
  sufficient decrease*" in ACM Transactions on Mathematical Software (TOMS)
  Volume **20**, Issue 3, Pages 286-307 (September 1994).

* Éric Thiébaut, "*Optimization issues in blind deconvolution algorithms*",
  SPIE Conf. Astronomical Data Analysis II, **4847**, 174-183 (2002).


## Installation

To install the library (not needed for Yorick nor for IDL):

1. Edit the file `src/Makefile` (see "Portability Issues" below).

2. Compile the library:

        cd src
        make

3. Optionaly install the software, for instance:

        make PREFIX=/usr/local install

   which will copy:

        liboptimpacklegacy.a into /usr/local/lib
        optimpacklegacy.h    into /usr/local/include

   and creates directory `/usr/local/doc/OptimPack-${VERSION}` with some
   documentation and legal stuff.


## Yorick Installation

### Building in any directory

Using the `configure` script in the `yorick` directory of OptimPackLegacy, it
is possible to build the Yorick plug-in in any directory (not necessarily the
source directory of the plug-in).

1. Make a new directory, move into that directory and run the `configure`
   script:

        mkdir -p $BUILD
        cd $BUILD
        $SRCDIR/configure

   where `$BUILD` is the name of the directory where to build the plug-in
   while`$SRCDIR` is the path to `yorick` directory of OptimPackLegacy.
   The `configure` script has a number of options do:

        $SRCDIR/configure --help

   for a description of the available options.


2. Compile the plugin code:

        make clean
        make

3. Optionaly install plugin in Yorick tree:

        make install


### Building in the source directory (not recommended)

1. Go to directory `yorick` of OptimPackLegacy.

2. Setup for compilation and compile the plugin code:

        yorick -batch make.i
        make clean
        make

3. Optionaly install plugin in Yorick tree:

        make install


## Portability Issues

OptimPackLegacy is written in standard ANSI-C99 and should pose no problems of
portability.
