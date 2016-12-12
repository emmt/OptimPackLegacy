/*
 * OptimPack1.i --
 *
 * Main startup file for OptimPack extension of Yorick.
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (c) 2003-2009, 2016 Éric Thiébaut.
 *
 * This file is part of OptimPack <https://github.com/emmt/OptimPackLegacy>.
 *
 * OptimPack is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option)
 * any later version.
 *
 * OptimPack is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * OptimPack (file "LICENSE" in the top source directory); if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307 USA
 *
 *-----------------------------------------------------------------------------
 */

if (is_func(plug_in)) plug_in, "OptimPack1";

extern __op_csrch;
/* PROTOTYPE
   int op_csrch(double f, double g, double array stp,
                double ftol, double gtol, double xtol,
		double stpmin, double stpmax, int array task,
		char array csave, long array isave, double array dsave);
*/

func op_csrch(f, g, &stp, ftol, gtol, xtol, stpmin, stpmax, &task,
              &csave, isave, dsave)
{
  if (numberof(task) != 1) error, "TASK must be a scalar";
  if (structof(isave) != long || numberof(isave) < 2)
    error, "bad ISAVE array";
  if (structof(dsave) != double || numberof(dsave) < 12)
    error, "bad DSAVE array";
  itask = int(task);
  cbuf = array(char, 128);
  info = __op_csrch(f, g, stp, ftol, gtol, xtol, stpmin, stpmax, itask,
                    cbuf, isave, dsave);
  task = long(itask);
  csave = string((task == 1 ? 0 : &cbuf));
  return long(info);
}

func op_vmlmb_setup(n, m, fmin=, fatol=, frtol=, sftol=, sgtol=, sxtol=,
                    epsilon=, costheta=)
{
  csave = array(char, 128);
  isave = array(long,  12);
  dsave = array(double, 27 + n + 2*m*(n + 1));
  if (is_void(frtol)) frtol = 1e-10;
  if (is_void(fatol)) fatol = 1e-13;
  if (is_void(sftol)) sftol = 0.001;
  if (is_void(sgtol)) sgtol = 0.9;
  if (is_void(sxtol)) sxtol = 0.1;
  if (is_void(epsilon)) epsilon = 1E-8;
  if (is_void(costheta)) costheta = 1E-2;
  task = long(__op_vmlmb_first(n, m, fatol, frtol, sftol, sgtol, sxtol,
                               epsilon, costheta, csave, isave, dsave));
  if (task != 1) error, string(&csave);
  ws = [&csave, &isave, &dsave];
  if (! is_void(fmin)) {
    op_vmlmb_set_fmin, ws, fmin;
  }
  return ws;
}

func op_vmlmb_msg(ws) { return string(ws(1)); }

func op_vmlmb_next(x, &f, &g, ws, active, h)
{
  local csave; eq_nocopy, csave, *ws(1);
  if (structof(csave) != char || numberof(csave) != 128)
    error, "corrupted workspace (CSAVE)";

  local isave; eq_nocopy, isave, *ws(2);
  if (structof(isave) != long || numberof(isave) != 12)
    error, "corrupted workspace (ISAVE)";
  m = isave(5);
  n = isave(6);

  local dsave; eq_nocopy, dsave, *ws(3);
  if (structof(dsave) != double || numberof(dsave) != 27 + n + 2*m*(1 + n))
    error, "corrupted workspace (DSAVE)";

  if (structof(x) != double || numberof(x) != n)
    error, "bad parameter array X";
  if (structof(f) != double || dimsof(f)(1))
    error, "bad function value F";
  if (structof(g) != double || numberof(g) != n)
    error, "bad gradient array G";

  if (! is_void(active) && (structof(active) != int ||
                            numberof(active) != n )) error, "bad array ACTIVE";

  if (! is_void(h) && (structof(h) != double ||
                       numberof(h) != n )) error, "bad array H";

  return long(__op_vmlmb_next(x, f, g, &active, &h, csave, isave, dsave));
}

extern __op_vmlmb_next;
/* PROTOTYPE
   int op_vmlmb_next(double array x, double array f, double array g,
                     pointer active, pointer h,
		     char array csave, long array isave, double array dsave);
*/

//extern op_vmlmb_first;
///* DOCUMENT op_vmlmb_first(n, m, fatol, frtol, sftol, sgtol, sxtol,
//                           csave, isave, dsave);
//     The returned value should be 1 unless there is an error.
//*/

extern __op_vmlmb_first;
/* PROTOTYPE
   int op_vmlmb_first(long n, long m,
                      double fatol, double frtol,
                      double sftol, double sgtol, double sxtol,
                      double epsilon, double costheta,
                      char array csave, long array isave, double array dsave);
*/

extern __op_vmlmb_set_fmin;
/* PROTOTYPE
   int op_vmlmb_set_fmin(char array csave,
                         long array isave,
                         double array dsave,
                         double new_value,
                         double array old_value);
*/
extern __op_vmlmb_get_fmin;
/* PROTOTYPE
   int op_vmlmb_get_fmin(char array csave,
                         long array isave,
                         double array dsave,
                         double array value);
*/

local op_vmlmb_set_fmin;
local op_vmlmb_get_fmin;
/* DOCUMENT old = op_vmlmb_set_fmin(ws, fmin);
         or fmin = op_vmlmb_get_fmin(ws);

      The function op_vmlmb_set_fmin set the value of FMIN in workspace WS and
      returns the previous value of FMIN if any (nil otherwise).

      The function op_vmlmb_get_fmin returns the actual value of FMIN in
      workspace WS or nil if FMIN has never been set.

   SEE ALSO: op_vmlmb_first.
 */
func op_vmlmb_set_fmin(ws, new)
{
  old = 0.0;
  if (__op_vmlmb_set_fmin(*ws(1), *ws(2), *ws(3), new, old)) {
    return old;
  }
}
func op_vmlmb_get_fmin(ws)
{
  fmin = 0.0;
  if (__op_vmlmb_get_fmin(*ws(1), *ws(2), *ws(3), fmin)) {
    return fmin;
  }
}

func op_mnb(f, x, &fx, &gx, fmin=, extra=,
            xmin=, xmax=, method=, mem=, verb=, quiet=,
            viewer=, printer=,
            maxiter=, maxeval=, output=,
            frtol=, fatol=, sftol=, sgtol=, sxtol=)
/* DOCUMENT op_mnb(f, x);
         or op_mnb(f, x, fout, gout);

     This is a wrapper to preserve compatibility with previous versions of
     OptimPack. See `op_vmlmb` for a full documentation.

   SEE ALSO op_vmlmb.
*/
{
  /* Choose minimization method.
     FIXME: METHOD - Scalar integer which  defines the optimization method to use.
     FIXME:     Conjugate  gradient   algorithm  is  used  if  one   of  the  bits
     FIXME:     OP_FLAG_POLAK_RIBIERE,         OP_FLAG_FLETCHER_REEVES,         or
     FIXME:     OP_FLAG_HESTENES_STIEFEL  is  set;  otherwise,  a  limited  memory
     FIXME:     variable  metric algorithm  (VMLM-B) is  used.  If  METHOD  is not
     FIXME:     specified and  if MEM=0, a conjugate gradient  search is attempted
     FIXME:     with flags: (OP_FLAG_UPDATE_WITH_GP |
     FIXME:                  OP_FLAG_SHANNO_PHUA    |
     FIXME:                  OP_FLAG_MORE_THUENTE   |
     FIXME:                  OP_FLAG_POLAK_RIBIERE  |
     FIXME:                  OP_FLAG_POWELL_RESTART)
     FIXME:     otherwise VMLM-B is used with flags: (OP_FLAG_UPDATE_WITH_GP |
     FIXME:                                           OP_FLAG_SHANNO_PHUA    |
     FIXME:                                           OP_FLAG_MORE_THUENTE).
     FIXME:     See documentation  of op_get_flags to  figure out the  allowed bit
     FIXME:     flags and their meaning.
  */
  if (! method) {
    /* Variable metric. */
    if (is_void(mem)) mem = min(n, 7);
    method = 0;
    method_name = swrite(format="Limited Memory BFGS (VMLM with MEM=%d)",
                         mem);
    ws = op_vmlmb_setup(n, mem, fmin=fmin,
                        fatol=fatol, frtol=frtol,
                        sftol=sftol, sgtol=sgtol, sxtol=sxtol);
  } else if (method < 0) {
    if (is_void(mem)) mem = min(n, 7);
    method_name = swrite(format="Limited Memory BFGS (LBFGS with MEM=%d)",
                         mem);
    ws = op_lbfgs_setup(n, mem);
  } else if (method >= 1 && method <= 15) {
    /* Conjugate gradient. */
    mem = 2;
    error, "conjugate-gradient method not yet implemented";
    method_name = swrite(format="Conjugate Gradient (%s)",
                         ["Fletcher-Reeves", "Polak-Ribiere",
                          "Polak-Ribiere with non-negative BETA"](method&3));
    ws = optim_cgmn_setup(method, fmin=fmin, fatol=fatol, frtol=frtol);
  } else {
    error, "bad METHOD";
  }
  return op_vmlmb(f, x, fx, gx, fmin=fmin, extra=extra,
                  xmin=xmin, xmax=xmax, flags=method, mem=mem,
                  verb=verb, quiet=quiet, viewer=viewer, printer=printer,
                  maxiter=maxiter, maxeval=maxeval, output=output,
                  frtol=frtol, fatol=fatol,
                  sftol=sftol, sgtol=sgtol, sxtol=sxtol);
}

func op_vmlmb(f, x, &fx, &gx, fmin=, extra=, xmin=, xmax=, flags=, mem=,
              verb=, quiet=, viewer=, printer=, maxiter=, maxeval=, output=,
              frtol=, fatol=, gatol=, grtol=, sftol=, sgtol=, sxtol=, savebest=)
/* DOCUMENT op_mnb(f, x);
         or op_mnb(f, x, fout, gout);

     Returns a minimum of a multivariate function by an iterative minimization
     algorithm (limited memory variable metric) possibly with simple bound
     constraints on the parameters.  Arguments are:

       F - User defined function to optimize.
           The prototype of F is:
             func F(x, &gx) {
               fx = ....; // compute function value at X
               gx = ....; // store gradient of F in GX
               return fx; // return F(X)
             }

       X - Starting solution (a floating point array).

       FOUT - Optional output variable to store the value of F at the
           minimum.

       GOUT - optional output variable to store the value of the gradient
           of F at the minimum.

     If the multivariate function has more than one minimum, which minimum
     is returned is undefined (although it depends on the starting
     parameters X).

     In case of early termination, the best solution found so far is
     returned.


   KEYWORDS

     EXTRA - Supplemental argument for F; if non-nil, F is called as
         F(X,GX,EXTRA) so its prototype must be: func F(x, &gx, extra).

     XMIN, XMAX  - Lower/upper bounds for  X.  Must be  conformable with X.
         For instance with XMIN=0, the non-negative solution will be
         returned.

     MEM - Number of previous directions used in variable metric limited
         memory method (default min(7, numberof(X))).

     MAXITER - Maximum number of iterations (default: no limits).

     MAXEVAL - Maximum number of function evaluations (default: no limits).

     FATOL, FRTOL - Relative function change tolerance for convergence (default:
         1.5e-8).

     GATOL, GRTOL - Absolute and relative gradient tolerances for convergence
         which is assumed whenever the Euclidean norm of the (projected)
         gradient is smaller or equal max(GATOL, GRTOL*GINIT) where GINIT is
         the Euclidean norm of the (projected) initila gradient.  By default,
         GTAOL=0 and GRTOL=1e-6.

     SAVEBEST - The default is to save the best variables among all tried ones
         in order to always return the best solution encountered by the
         algorithm.  This has however a memory cost, you may set this option
         explicitly false (zero) to save a bit of memory, but in case of early
         stopping, the returned solution may not be the best one.

     VERB - Verbose mode?  If non-nil and non-zero, print out information
         every VERB iterations and for the final one.

     QUIET - If true and not in verbose mode, do not print warning nor
         convergence error messages.

     OUPTPUT - Output for verbose mode.  For instance, text file stream
         opened for writing.

     VIEWER - User defined subroutine to call every VERB iterations (see
         keyword VERB above)to display the solution X.  The subroutine will be
         called as:

            viewer, x, extra;

         where X is the current solution and EXTRA is the value of keyword
         EXTRA (which to see).  If the viewer uses Yorick graphics window(s)
         it may call "pause, 1;" before returning to make sure that graphics
         get correctly updated.

     PRINTER - User defined subroutine to call every VERB iterations (see
         keyword VERB above) to printout iteration information.  The
         subroutine will be called as:

            printer, output, iter, eval, cpu, fx, gnorm, steplen, x, extra;

         where OUTPUT is the value of keyword OUTPUT (which to see), ITER is
         the number of iterations, EVAL is the number of function evaluations,
         CPU is the elapsed CPU time in seconds, FX is the function value at
         X, GNORM is the Euclidean norm of the gradient at X, STEPLEN is the
         length of the step along the search direction, X is the current
         solution and EXTRA is the value of keyword EXTRA (which to see).

     SFTOL, SGTOL, SXTOL - Line search tolerance and safeguard
        parameters (see op_csrch).

   SEE ALSO: op_get_flags, op_csrch,
             op_vmlmb_setup, op_vmlmb_next.
 */
{
  local result, gx;

  /* Get function. */
  if (! is_func(f)) {
    error, "expecting a function for argument F";
  }
  use_extra = (! is_void(extra));

  /* Starting parameters. */
  if ((s = structof(x)) != double && s != float && s != long &&
      s != int && s != short && s != char) {
    error, "expecting a numerical array for initial parameters X";
  }
  n = numberof(x);
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
  check_iter = (! is_void(maxiter));
  check_eval = (! is_void(maxeval));
  if (check_eval && maxeval < 1) {
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

  /* By default, save the best solution so far. */
  if (is_void(savebest)) savebest = 1n;

  /* Choose minimization method. */
  if (is_void(mem)) mem = min(n, 7);
  method_name = swrite(format="VMLMB %s bounds and MEM=%d",
                       (bounds != 0 ? "with" : "without"), mem);
  ws = op_vmlmb_setup(n, mem, fmin=fmin,
                      fatol=fatol, frtol=frtol,
                      sftol=sftol, sgtol=sgtol, sxtol=sxtol);

  /* Start iterations. */
  step = 0.0;
  task = 1;
  eval = iter = 0;
  stop = 0n;
  if (verb) {
    elapsed = array(double, 3);
    timer, elapsed;
    cpu_start = elapsed(1);
  }
  if (structof(x) != double) {
    x = double(x);
  }
  local best_x, best_fx, best_gnorm; /* best solution so far */
  local gx, gnorm; /* the gradient and its (projected) norm */
  local active;
  for (;;) {
    if (task == 1) {
      /* Evaluate function and gradient. */
      if (check_eval && eval >= maxeval) {
        stop = 1n;
        msg = swrite(format="warning: too many function evaluations (%d)\n",
                     eval);
      } else {
        if (bounds != 0) {
          if ((bounds & 1) == 1) {
            x = max(x, xmin);
          }
          if ((bounds & 2) == 2) {
            x = min(x, xmax);
          }
        }
        fx = (use_extra ? f(x, gx, extra) : f(x, gx));
        ++eval;
        if (bounds != 0) {
          /* Figure out the set of free parameters:
           *   ACTIVE(i) = 0 if X(i) has a lower bound XMIN(i)
           *                 and X(i) = XMIN(i) and GX(i) >= 0
           *               0 if X(i) value has an upper bound XMAX(i)
           *                 and X(i) = XMAX(i) and GX(i) <= 0
           *               1 (or any non-zero value) otherwise
           */
          if (bounds == 1) {
            active = ((x > xmin) | (gx < 0.0));
          } else if (bounds == 2) {
            active = ((x < xmax) | (gx > 0.0));
          } else {
            active = (((x > xmin) | (gx < 0.0)) & ((x < xmax) | (gx > 0.0)));
          }
        }
        if (eval == 1 && grtol > 0) {
          gnorm = (bounds != 0 ? op_norm2(active*gx) : op_norm2(gx));
          gtest = max(gtest, grtol*gnorm);
        } else {
          gnorm = -1; // to indicate that it must be updated if needed
        }
        if (savebest && (eval == 1 || fx < best_fx)) {
          best_x = x;
          best_fx = fx;
          if (gnorm < 0) {
            gnorm = (bounds != 0 ? op_norm2(active*gx) : op_norm2(gx));
          }
          best_gnorm = gnorm;
        }
      }
    }

    /* Check for convergence. */
    if (task != 1 || eval == 1) {
      if (gnorm < 0) {
        gnorm = (bounds != 0 ? op_norm2(active*gx) : op_norm2(gx));
      }
      if (task > 2) {
        stop = 1n;
        msg = op_vmlmb_msg(ws);
      } else if (gnorm <= gtest) {
        stop = 1n;
        msg = swrite(format="convergence (%s)\n", "gradient small enough");
      } else if (check_iter && iter > maxiter) {
        stop = 1n;
        msg = swrite(format="warning: too many iterations (%d)\n", iter);
      }
      if (stop && savebest) {
        /* restore best solution so far. */
        eq_nocopy, x, best_x;
        fx = best_fx;
        gnorm = best_gnorm;
      }
      if (verb) {
        if (eval == 1 && ! use_printer) {
          write, output, format="# Method %d (MEM=%d): %s\n#\n",
            method, mem, method_name;
          write, output, format="# %s\n# %s\n",
            "ITER  EVAL   CPU (ms)        FUNC               GNORM   STEPLEN",
            "---------------------------------------------------------------";
        }
        if (stop || ! (iter % verb)) {
          timer, elapsed;
          cpu = 1e3*(elapsed(1) - cpu_start);
          if (use_printer) {
            printer, output, iter, eval, cpu, fx, gnorm, steplen, x, extra;
          } else {
            write, output, format=" %5d %5d %10.3f  %+-24.15e%-9.1e%-9.1e\n",
              iter, eval, cpu, fx, gnorm, step;
          }
          if (use_viewer) {
            viewer, x, extra;
          }
        }
      }
      if (stop) {
        if (msg && (verb || (task != 3 && ! quiet))) {
          write, output, format="# %s\n", strtrim(msg, 2, blank=" \t\v\n\r");
        }
        return x;
      }
    }

    /* Call optimizer. */
    task = op_vmlmb_next(x, fx, gx, ws, active);
    iter = (*ws(2))(7);
    step = (*ws(3))(22);
  }
}

func op_norm2(x) { return sqrt(sum(x*x)); }
/* DOCUMENT op_norm2(x);

     Yield the Euclidean norm of X.
*/
