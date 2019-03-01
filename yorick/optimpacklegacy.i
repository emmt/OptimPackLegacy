/*
 * optimpacklegacy.i --
 *
 * Main startup file for OptimPackLegacy extension of Yorick.
 *
 *-----------------------------------------------------------------------------
 *
 * This file is part of OptimPackLegacy
 * <https://github.com/emmt/OptimPackLegacy>.
 *
 * Copyright (c) 2003-2009, 2016 Éric Thiébaut.
 *
 * OptimPackLegacy is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option)
 * any later version.
 *
 * OptimPackLegacy is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * OptimPackLegacy (file "LICENSE" in the top source directory); if not, write
 * to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA 02111-1307 USA
 *
 *-----------------------------------------------------------------------------
 */

if (is_func(plug_in)) plug_in, "optimpacklegacy";

local OPL_TASK_START, OPL_TASK_FG, OPL_TASK_FREEVARS, OPL_TASK_NEWX;
local OPL_TASK_CONV, OPL_TASK_WARN, OPL_TASK_ERROR;
extern opl_vmlmb_create;
/* DOCUMENT ws = opl_vmlmb_create(dims, mem, key1=val1, key2=val2, ...);

     Create a new workspace for the VMLMB algorithm.  DIMS gives the dimension
     list of the variables (like `dimsof`) and MEM is the number of previous
     steps to memorize.

     The workspace WS can be used to retrieve the following attributes:

       ws.dims ............. Dimension list of the variables.
       ws.size ............. Size of the problem (number of variables).
       ws.mem .............. Number of previous steps to memorize.
       ws.task ............. Current pending task.
       ws.evaluations ...... Number of function calls.
       ws.iterations ....... Number of iterations.
       ws.restarts ......... Number of algorithm restarts.
       ws.step ............. Lenght of the current or last step.
       ws.status ........... Current status value.
       ws.reason ........... Explanatory message about current status/task.
       ws.gnorm ............ Euclidean norm of the (projected) gradient at the
                             last successful step.
       ws.fmin ............. Strict lower bound for the function (NaN if not
                             set).
       ws.fatol ............ Absolute function tolerance for the convergence.
       ws.frtol ............ Relative function tolerance for the convergence.
       ws.sftol ............ Tolerance for the sufficient decrease condition.
       ws.sgtol ............ Tolerance for the curvature condition.
       ws.sxtol ............ Relative tolerance for an acceptable step.
       ws.delta ............ Relative size of a small step.
       ws.epsilon .......... Threshold for the sufficient descent condition.

     Attributes FMIN, FATOL, FRTOL, SFTOL, SGTOL, SXTOL, DELTA and EPSILON are
     configurable via keywords when calling `opl_vmlmb_create`, they may be
     configure later with `opl_vmlmb_configure`.  The other attributes are
     read-only.

     It not recommended to have MEM larger than the number of variables and it
     is often the case that a modest value for MEM, say MEM=5, is as efficient
     as larger values.  As you may expect, MEM has an incidence on the size of
     the neessary memory which is about `(2*m*(n + 1) + n)*sizeof(double)`
     where `n` is the number of variables and `m` is MEM.  Allocated memory is
     automatically released when the returned workspace is no longer in use.

     The pending task can take one of the following values:

       OPL_TASK_START ...... Start line search.
       OPL_TASK_FG ......... Caller has to compute function and gradient.
       OPL_TASK_FREEVARS ... Caller has to determine the free variables.
       OPL_TASK_NEWX ....... New variables available for inspection.
       OPL_TASK_CONV ....... Search has converged.
       OPL_TASK_WARN ....... Search aborted with warning.
       OPL_TASK_ERROR ...... Search aborted with error.


   SEE ALSO opl_vmlmb_configure, opl_vmlmb_iterate, dimsof.
*/

extern opl_vmlmb_configure;
/* DOCUMENT opl_vmlmb_configure, ws, key1=val1, key2=val2, ...;

     Configure VMLMB workspace WS.  Parameters are provided as keyword-value
     pairs.  See `opl_vmlmb_create` for a list of configurable settings. When
     called as a function, WS is returned.

   SEE ALSO opl_vmlmb_create
*/

extern opl_vmlmb_iterate;
/* DOCUMENT task = opl_vmlmb_iterate(ws, x, f, g);
         or task = opl_vmlmb_iterate(ws, x, f, g, isfree);
         or task = opl_vmlmb_iterate(ws, x, f, g, isfree, h);

     Perform one step of the VMLMB algorithm.  WS is VMLMB workspace, X gives
     the variables, F and G are the function value and gradient at X.  ISFREE
     is an optional array which indicates which variables are free to vary.
     H is an optional array which provides a diagonal preconditioner.  The
     returned value is the next pending task.

     All specified arrays must have the same dimensions as those expected by
     WS.  F must be a simple variable reference (not an expression), X and G
     must be arrays of double's; if specified, ISFREE is an array of int's; if
     specified, H is an array of double's.

   SEE ALSO opl_vmlmb_create, opl_vmlmb_restore.
*/

extern opl_vmlmb_restore;
/* DOCUMENT task = opl_vmlmb_restore(ws, x, f, g);

     Restore last line search starting point for VMLMB workspace WS.  Calling
     this is only effective if task is OPL_TASK_FG.  Arguments X, F and G are
     the same as in `opl_vmlmb_iterate`.

   SEE ALSO opl_vmlmb_iterate.
*/

extern opl_vmlmb_restart;
/* DOCUMENT task = opl_vmlmb_restart(ws);

      Set VMLMB workspace WS so that it can be used for a new optimization with
      the same parameters.

   SEE ALSO opl_vmlmb_create, opl_vmlmb_iterate.
*/

func opl_vmlmb(f, x, &fx, &gx, fmin=, extra=, xmin=, xmax=, flags=, mem=,
               verb=, quiet=, viewer=, printer=, maxiter=, maxeval=, output=,
               frtol=, fatol=, gatol=, grtol=, sftol=, sgtol=, sxtol=)
/* DOCUMENT opl_vmlmb(f, x);
         or opl_vmlmb(f, x, fout, gout);

     Returns a minimum of a multivariate function F by an iterative
     minimization algorithm (limited memory variable metric) possibly with
     simple bound constraints on the parameters.  F is the function
     to minimize, its prototype is:

         func f(x, &gx) {
             fx = ....; // compute function value at X
             gx = ....; // store gradient of F in GX
             return fx; // return F(X)
         }

     Argument X is the starting solution (a double precision floating point
     array).  FOUT and GOUT are optional output variables to store the value of
     F and its gradient at the minimum.

     If the multivariate function has more than one minimum, which minimum is
     returned is undefined (although it depends on the starting parameters X).

     In case of early termination, the best solution found so far is returned.


   KEYWORDS

     EXTRA - Supplemental argument for F; if non-nil, F is called as
         F(X,GX,EXTRA) so its prototype must be: func F(x, &gx, extra).  It is
         however more flexible to use a closure for F if additional data is
         needed.

     XMIN, XMAX - Lower/upper bounds for X.  Must be conformable with X.  For
         instance with XMIN=0, the non-negative solution will be returned.

     MEM - Number of previous directions used in variable metric limited memory
         method (default min(7, numberof(X))).

     MAXITER - Maximum number of iterations (default: no limits).

     MAXEVAL - Maximum number of function evaluations (default: no limits).

     FATOL, FRTOL - Relative function change tolerance for convergence
         (default: 1.5e-8).

     GATOL, GRTOL - Absolute and relative gradient tolerances for convergence
         which is assumed whenever the Euclidean norm of the (projected)
         gradient is smaller or equal max(GATOL, GRTOL*GINIT) where GINIT is
         the Euclidean norm of the (projected) initial gradient.  By default,
         GATOL=0 and GRTOL=1e-6.

     VERB - Verbose mode?  If non-nil and non-zero, print out information every
         VERB iterations and for the final one.

     QUIET - If true and not in verbose mode, do not print warning nor
         convergence error messages.

     OUPTPUT - Output for verbose mode.  For instance, text file stream opened
         for writing.

     VIEWER - User defined subroutine to call every VERB iterations (see
         keyword VERB above)to display the solution X.  The subroutine will be
         called as:

            viewer, x, extra;

         where X is the current solution and EXTRA is the value of keyword
         EXTRA (which to see).  If the viewer uses Yorick graphics window(s) it
         may call "pause, 1;" before returning to make sure that graphics get
         correctly updated.

     PRINTER - User defined subroutine to call every VERB iterations (see
         keyword VERB above) to printout iteration information.  The subroutine
         will be called as:

            printer, output, iter, eval, cpu, fx, gnorm, steplen, x, extra;

         where OUTPUT is the value of keyword OUTPUT (which to see), ITER is
         the number of iterations, EVAL is the number of function evaluations,
         CPU is the elapsed CPU time in seconds, FX is the function value at X,
         GNORM is the Euclidean norm of the gradient at X, STEPLEN is the
         length of the step along the search direction, X is the current
         solution and EXTRA is the value of keyword EXTRA (which to see).

     SFTOL, SGTOL, SXTOL - Line search tolerance and safeguard parameters (see
        opl_csrch).

   SEE ALSO: opl_vmlmb_config, opl_vmlmb_create, opl_vmlmb_iterate.
*/
{
  /* Largest value of a long integer. */
  LONG_MAX = (1 << (sizeof(long)*8 - 1)) - 1;

  /* Get function. */
  if (is_void(f) || is_array(f)) {
    error, "expecting a function for argument F";
  }
  use_extra = (! is_void(extra));

  /* Starting parameters. */
  if ((s = structof(x)) != double && s != float && s != long &&
      s != int && s != short && s != char) {
    error, "expecting a numerical array for initial parameters X";
  }
  dims = dimsof(x);

  /* Bounds on parameters. */
  bounds = 0;
  if (! is_void(xmin)) {
    if (is_void((t = dimsof(x, xmin))) || t(1) != dims(1)
        || anyof(t != dims)) {
      error, "bad dimensions for lower bound XMIN";
    }
    if ((convert = (s = structof(xmin)) != double) && s != float &&
        s != long && s != int && s != short && s != char) {
      error, "bad data type for lower bound XMIN";
    }
    if (convert || (t = dimsof(xmin))(1) != dims(1) || anyof(t != dims)) {
      xmin += array(double, dims);
    }
    bounds |= 1;
  }
  if (! is_void(xmax)) {
    if (is_void((t = dimsof(x, xmax))) || t(1) != dims(1)
        || anyof(t != dims)) {
      error, "bad dimensions for lower bound XMAX";
    }
    if ((convert = (s = structof(xmax)) != double) && s != float &&
        s != long && s != int && s != short && s != char) {
      error, "bad data type for lower bound XMAX";
    }
    if (convert || (t = dimsof(xmax))(1) != dims(1) || anyof(t != dims)) {
      xmax += array(double, dims);
    }
    bounds |= 2;
  }

  /* Output stream. */
  if (! is_void(output)) {
    if (structof(output) == string) {
      output = open(output, "a");
    } else if (typeof(output) != "text_stream") {
      error, "bad value for keyword OUTPUT";
    }
  }

  /* Maximum number of iterations and function evaluations. */
  if (is_void(maxiter)) maxiter = LONG_MAX;
  if (is_void(maxeval)) maxeval = LONG_MAX;
  if (maxeval < 1) {
    error, "MAXEVAL must be at least 1";
  }

  /* Viewer and printer subroutines. */
  use_printer = (! is_void(printer));
  use_viewer  = (! is_void(viewer));

  /* Global convergence parameters. */
  if (is_void(gatol)) {
    gatol = 0.0;
  } else if (is_scalar(gatol) && identof(gatol) <= Y_DOUBLE &&
             gatol >= 0) {
    gatol = double(gatol);
  } else {
    error, "bad value for GATOL";
  }
  if (is_void(grtol)) {
    grtol = 0.0;
  } else if (is_scalar(grtol) && identof(grtol) <= Y_DOUBLE &&
             grtol >= 0 && grtol <= 1) {
    grtol = double(grtol);
  } else {
    error, "bad value for GRTOL";
  }
  gtest = double(gatol);

  /* Choose minimization method. */
  if (is_void(mem)) mem = min(numberof(x), 7);
  method_name = swrite(format="VMLMB %s bounds and MEM=%d",
                       (bounds != 0 ? "with" : "without"), mem);
  ws = opl_vmlmb_create(dims, mem, fmin=fmin,
                        fatol=fatol, frtol=frtol,
                        sftol=sftol, sgtol=sgtol, sxtol=sxtol);

  /* Start iterations. */
  task = 1;
  eval = 0;
  stop = 0n;
  if (verb) {
    elapsed = array(double, 3);
    timer, elapsed;
    cpu_start = elapsed(1);
    print_header = (! use_printer);
    last_print_iter = -1;
  }
  if (structof(x) != double) {
    x = double(x);
  }
  local gx, gnorm, isfree, iter, step;
  task = ws.task;
  for (;;) {
    if (task == OPL_TASK_FG) {
      /* Evaluate function and gradient. */
      if (eval >= maxeval) {
        /* Too many function evaluations.  We restore the variables at the
           start of the line search which is a cheap way (no extra memory cost)
           to recover variables which should be nearly the best ones.  We wait
           until last iteration information have been printed to stop the
           iterations. */
        stop = 1n;
        msg = swrite(format="too many function evaluations (%d)", eval);
        opl_vmlmb_restore, ws, x, fx, gx;
      } else {
        if (bounds != 0) {
          if ((bounds & 1) == 1) {
            x = max(unref(x), xmin);
          }
          if ((bounds & 2) == 2) {
            x = min(unref(x), xmax);
          }
        }
        fx = (use_extra ? f(x, gx, extra) : f(x, gx));
        ++eval;
      }
    }
    if (task == OPL_TASK_FREEVARS && bounds != 0) {
      /* Determine the set of free variables. */
      isfree = [];
      if (bounds == 1) {
        isfree = ((x > xmin) | (gx < 0.0));
      } else if (bounds == 2) {
        isfree = ((x < xmax) | (gx > 0.0));
      } else {
        isfree = (((x > xmin) | (gx < 0.0)) & ((x < xmax) | (gx > 0.0)));
      }
    }

    /* Check for convergence. */
    if (task >= OPL_TASK_NEWX) {
      iter = ws.iterations;
      if (task >= OPL_TASK_WARN) {
        /* Error or warning. */
        stop = 1n;
        msg = ws.reason;
      } else if (ws.gnorm <= gtest) {
        stop = 1n;
        msg = swrite(format="convergence (%s)", "gradient small enough");
      } else if (iter >= maxiter) {
        stop = 1n;
        msg = swrite(format="too many iterations (%d)", iter);
      }
    }
    if (verb && (stop || task >= OPL_TASK_NEWX && (iter % verb) == 0) &&
        iter > last_print_iter) {
      /* Print information. */
      if (print_header) {
        write, output, format="# Method %d (MEM=%d): %s\n#\n",
          0, mem, method_name;
        write, output, format="# %s\n# %s\n",
          "ITER  EVAL   CPU (ms)        FUNC               GNORM   STEPLEN",
          "---------------------------------------------------------------";
        print_header = 0n;
      }
      timer, elapsed;
      cpu = 1e3*(elapsed(1) - cpu_start);
      step = ws.step;
      gnorm = ws.gnorm;
      if (use_printer) {
        printer, output, iter, eval, cpu, fx, gnorm, step, x, extra;
      } else {
        write, output, format=" %5d %5d %10.3f  %+-24.15e%-9.1e%-9.1e\n",
          iter, eval, cpu, fx, gnorm, step;
      }
      if (use_viewer) {
        viewer, x, extra;
      }
      last_print_iter = iter;
    }
    if (stop) {
      if (msg && (verb || (task != 3 && ! quiet))) {
        write, output, format="# %s\n", strtrim(msg, 2, blank=" \t\v\n\r");
      }
      return x;
    }

    /* Call optimizer. */
    task = opl_vmlmb_iterate(ws, x, fx, gx, isfree);
  }
}

extern _opl_init;
/* DOCUMENT _opl_init;

     Restore/set global variables used by YOPL.  In principle, it is not needed
     to call this subroutine (this is automatically done at startup) unless you
     destroyed some constants.
*/
_opl_init; /* restore global variables */
