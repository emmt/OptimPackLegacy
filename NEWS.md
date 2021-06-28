## Version 1.5.0 (2019-XX-XX)

* A Python interface to VMLM-B is provided (thanks to Romain Fetick).

* Array dimensions and indices are now `long`.

* Compilation and installation is much easier thanks to GNU-autoconf tools.

* With `VMLMB`, **warm starting** is possible to re-use the collected Hessian
  information when minimizing a (slightly) different objective function while
  the variables and the algorithm settings remain the same.

* Scaling parameter `gamma` is computed once if there are no constraints.

* Code is in standard C99.


## Version 1.4.0 (2017-01-10)

* VMLMB workspace is now a single opaque handle.  All needed data can be saved
  in a single block of memory (*monolithic* workspace) provided by the caller
  or allocated by the library (and possibly split in small blocks).  This
  simplifies the use of VMLMB both in C and in other programming languages
  and add flexibility for very large problems.

* In case of early stopping (i.e., during a line search), it is now posible to
  restore the variables at the start of the line search.  This is expected to
  be nearly the best ones so far.

* Rewrite VMLMB code for better readability, to remove unecessary variables and
  avoid unecessary computations.  This also fixes some inconsistencies.  Opaque
  structures are now used instead of arrays of values, hopefully not to the
  detriment of portability.

* IDL port has been temporarily removed (will be back soon).

* All functions, constants and macros prefixes become `opl_` or `OPL_`.


## Version 1.3.3 (2016-12-13)

* Parameter `delta > 0` is used to specify a relative small step size.

* Parameter `epsilon ≥ 0` is used to check for the sufficient descent
  condition.  Seems to speed up convergence in many cases.

* In Yorick interface:

  - The norm of the *projected* gradient is displayed and used to check for the
    convergence when there are bounds.

  - New `op_vmlmb` interface in Yorick which controls the convergence based on
    the gradient via keywords `gatol` and `grtol`.

  - The `op_mnb` driver is deprecated in favor of `op_vmlmb`.

  - By default, the best solution found so far is returned by `op_vmlmb` (see
    keyword `savebest`).  convergence when there are bounds.

* Code cleanup: get rid of non-linear conjugate gradients (never implemented in
  this version).

* `make distrib VERSION=...` can be used to archive the sources.

* Move sources to `src` and update make target/rules.

* OptimPack (version 1) is now on GitHub at
  https://github.com/emmt/OptimPackLegacy


## Version 1.3.2 (2010-05-07)

* Makefile fixed thanks to Michael Williamson.


## Version 1.3.1 (2009-09-23)

* Bugs fixed in Yorick interface with the size of workspace arrays
  ISAVE and DSAVE.


## Version 1.3.0 (2009-05-14)

* `fmin` is now an optional parameter.  Functions `op_get/set_fmin` and
  `op_vmlmb_first` have changed.  Accordingly, the Yorick interface have
  changed.

* IDL interface needs to be updated...


## Version 1.2.0 (2008-01-31)


## Version 1.1.0 (2008-01-30)

* Yorick interface to OptimPack reworked.

* In IDL interface, some fixes by Laurent Mugnier:

  - syntax errors in `op_vmlmb_msg.pro`;

  - bugs in parsing of arguments in `op_idl.c` (function `op_idl_vmlmb_next`).

  - New function `fmin_op` in directory `idl/contrib` by Laurent Mugnier
    (license is CeCILL not GPL).

* GNU Public License updated.

* Integrate changes made by Laurent Mugnier (thanks to him!)  in IDL interface
  to support for Linux and 64 bit machines.


## Version 1.0.0 (2003-04-22)
