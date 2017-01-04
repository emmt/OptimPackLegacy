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

Directory [yorick](yorick) contains an implementation of OptimPack support in
Yorick.


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

1. Go to directory `yorick`.

2. Setup for compilation and compile the plugin code:

        yorick -batch make.i
        make clean
        make

3. Optionaly install plugin in Yorick tree:

        make install


## Portability Issues

OptimPackLegacy is written in standard ANSI-C and should pose no problem of
portability.  However, in order to match the data types used in your software,
you may have to set the values of the following macros:

    OPL_INTEGER = data type used to store array indices
    OPL_LOGICAL = data type of the result of a logical test

This must be  done *before* `optimpacklegacy.h` get included.   If these macros
are not defined, the default assumed in `optimpacklegacy.h` is:

    OPL_INTEGER = int
    OPL_LOGICAL = int

For instance, one should write:

    #define OPL_INTEGER long
    #define OPL_LOGICAL int
    #include "optimpacklegacy.h"
    ...

Of course, the installed OptimPackLegacy library must have been compiled with
the correct data types.  Another possibility is to define these macros when
calling CPP (the C preprocessor), *e.g.* in `src/Makefile`:

     CPPFLAGS = -DOPL_INTEGER=long -DOPL_LOGICAL=int -I.

a final possibility is to edit `optimpacklegacy.h` and to adjust the default
values of these macros (at the very beginning of this file).  If you plan to
install in your system, the best is probably to fix the definitions in
`optimpacklegacy.h`, then compile the library and finally install the library
and the header file `optimpacklegacy.h` in proper system directories (*e.g.*
with `make install PREFIX=...`).
