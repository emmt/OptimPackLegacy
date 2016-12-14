/*
 * optimpacklegacy.i --
 *
 * Main startup file for OptimPackLegacy extension of Yorick.
 *
 *-----------------------------------------------------------------------------
 *
 * This file is part of OptimPack <https://github.com/emmt/OptimPackLegacy>.
 *
 * Copyright (c) 2003-2009, 2016 Éric Thiébaut.
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

if (is_func(plug_in)) plug_in, "optimpacklegacy";

extern __opl_csrch;
/* PROTOTYPE
   int opl_csrch(double f, double g, double array stp,
                double ftol, double gtol, double xtol,
		double stpmin, double stpmax, int array task,
		char array csave, long array isave, double array dsave);
*/

func opl_csrch(f, g, &stp, ftol, gtol, xtol, stpmin, stpmax, &task,
               &csave, isave, dsave)
{
  if (numberof(task) != 1) error, "TASK must be a scalar";
  if (structof(isave) != long || numberof(isave) < 2)
    error, "bad ISAVE array";
  if (structof(dsave) != double || numberof(dsave) < 12)
    error, "bad DSAVE array";
  itask = int(task);
  cbuf = array(char, 128);
  info = __opl_csrch(f, g, stp, ftol, gtol, xtol, stpmin, stpmax, itask,
                    cbuf, isave, dsave);
  task = long(itask);
  csave = string((task == 1 ? 0 : &cbuf));
  return long(info);
}

func opl_vmlmb_setup(n, m, fmin=, fatol=, frtol=, sftol=, sgtol=, sxtol=,
                     delta=, epsilon=)
{
  csave = array(char, 128);
  isave = array(long,  11);
  dsave = array(double, 28 + n + 2*m*(n + 1));
  if (is_void(frtol)) frtol = 1e-10;
  if (is_void(fatol)) fatol = 1e-13;
  if (is_void(sftol)) sftol = 0.001; // FIXME:
  if (is_void(sgtol)) sgtol = 0.9;
  if (is_void(sxtol)) sxtol = 0.1;
  if (is_void(delta)) delta = 1e-3;
  if (is_void(epsilon)) epsilon = 0.0;
  task = long(__opl_vmlmb_setup(n, m, fatol, frtol, sftol, sgtol, sxtol,
                               delta, epsilon, csave, isave, dsave));
  if (task != 1) error, string(&csave);
  ws = [&csave, &isave, &dsave];
  if (! is_void(fmin)) {
    opl_vmlmb_set_fmin, ws, fmin;
  }
  return ws;
}

func opl_vmlmb_restart(ws)
{
  /* Extract workspace data. */
  local m, n, csave, isave, dsave;
  __opl_vmlmb_parse, ws;

  /* Call wrapper. */
  return long(__opl_vmlmb_restart(csave, isave, dsave));
}

func opl_vmlmb_restore(x, &f, g, ws)
{
  /* Extract workspace data. */
  local m, n, csave, isave, dsave;
  __opl_vmlmb_parse, ws;

  /* Check other arguments. */
  if (structof(x) != double || numberof(x) != n) {
    error, "bad parameter array X";
  }
  if (structof(g) != double || numberof(g) != n) {
    error, "bad gradient array G";
  }

  /* Call wrapper. */
  fp = [0.0];
  task = long(__opl_vmlmb_restore(x, fp, g, csave, isave, dsave));
  f = fp(1);
  return task;
}

func opl_vmlmb_next(x, &f, g, ws, isfree, h)
{
  /* Extract workspace data. */
  local m, n, csave, isave, dsave;
  __opl_vmlmb_parse, ws;

  /* Check other arguments. */
  if (structof(x) != double || numberof(x) != n) {
    error, "bad parameter array X";
  }
  if (identof(f) > Y_DOUBLE || ! is_scalar(f)) {
    error, "bad function value F";
  }
  if (structof(g) != double || numberof(g) != n) {
    error, "bad gradient array G";
  }
  if (! is_void(isfree) && (structof(isfree) != int ||
                            numberof(isfree) != n )) {
    error, "bad array ISFREE";
  }
  if (! is_void(h) && (structof(h) != double ||
                       numberof(h) != n )) {
    error, "bad array H";
  }

  /* Call wrapper. */
  fp = [double(f)]
  task = long(__opl_vmlmb_next(x, fp, g, &isfree, &h, csave, isave, dsave));
  f = fp(1);
  return task;
}

func opl_vmlmb_get_msg(ws)         { return string(ws(1)); }

func opl_vmlmb_get_task(ws)        { return (*ws(2))( 3); }
func opl_vmlmb_get_m(ws)           { return (*ws(2))( 4); }
func opl_vmlmb_get_n(ws)           { return (*ws(2))( 5); }
func opl_vmlmb_get_iterations(ws)  { return (*ws(2))( 6); }
func opl_vmlmb_get_mark(ws)        { return (*ws(2))( 7); }
func opl_vmlmb_get_mp(ws)          { return (*ws(2))( 8); }
func opl_vmlmb_get_flags(ws)       { return (*ws(2))( 9); }
func opl_vmlmb_get_evaluations(ws) { return (*ws(2))(10); }
func opl_vmlmb_get_restarts(ws)    { return (*ws(2))(11); }

func opl_vmlmb_get_sftol(ws)       { return (*ws(3))(13); }
func opl_vmlmb_get_sgtol(ws)       { return (*ws(3))(14); }
func opl_vmlmb_get_sxtol(ws)       { return (*ws(3))(15); }
func opl_vmlmb_get_frtol(ws)       { return (*ws(3))(16); }
func opl_vmlmb_get_fatol(ws)       { return (*ws(3))(17); }
func opl_vmlmb_get_fmin(ws)        { return (*ws(3))(18); }
func opl_vmlmb_get_f0(ws)          { return (*ws(3))(19); }
func opl_vmlmb_get_gd(ws)          { return (*ws(3))(20); }
func opl_vmlmb_get_gd0(ws)         { return (*ws(3))(21); }
func opl_vmlmb_get_step(ws)        { return (*ws(3))(22); }
func opl_vmlmb_get_stepmin(ws)     { return (*ws(3))(23); }
func opl_vmlmb_get_stepmax(ws)     { return (*ws(3))(24); }
func opl_vmlmb_get_delta(ws)       { return (*ws(3))(25); }
func opl_vmlmb_get_epsilon(ws)     { return (*ws(3))(26); }
func opl_vmlmb_get_gnorm(ws)       { return (*ws(3))(27); }
func opl_vmlmb_get_g0norm(ws)      { return (*ws(3))(28); }

local opl_vmlmb_set_fmin;
local opl_vmlmb_get_fmin;
/* DOCUMENT old = opl_vmlmb_set_fmin(ws, fmin);
         or fmin = opl_vmlmb_get_fmin(ws);

      The function opl_vmlmb_set_fmin set the value of FMIN in workspace WS and
      returns the previous value of FMIN if any (nil otherwise).

      The function opl_vmlmb_get_fmin returns the actual value of FMIN in
      workspace WS or nil if FMIN has never been set.

   SEE ALSO: opl_vmlmb_setup.
 */
func opl_vmlmb_set_fmin(ws, new)
{
  old = 0.0;
  if (__opl_vmlmb_set_fmin(*ws(1), *ws(2), *ws(3), new, old)) {
    return old;
  }
}
func opl_vmlmb_get_fmin(ws)
{
  fmin = 0.0;
  if (__opl_vmlmb_get_fmin(*ws(1), *ws(2), *ws(3), fmin)) {
    return fmin;
  }
}

func opl_vmlmb(f, x, &fx, &gx, fmin=, extra=, xmin=, xmax=, flags=, mem=,
              verb=, quiet=, viewer=, printer=, maxiter=, maxeval=, output=,
              frtol=, fatol=, gatol=, grtol=, sftol=, sgtol=, sxtol=)
/* DOCUMENT opl_vmlmb(f, x);
         or opl_vmlmb(f, x, fout, gout);

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
         the Euclidean norm of the (projected) initila gradient.  By default,
         GTAOL=0 and GRTOL=1e-6.

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

   SEE ALSO: opl_get_flags, opl_csrch,
             opl_vmlmb_setup, opl_vmlmb_next.
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
  if (is_void(mem)) mem = min(n, 7);
  method_name = swrite(format="VMLMB %s bounds and MEM=%d",
                       (bounds != 0 ? "with" : "without"), mem);
  ws = opl_vmlmb_setup(n, mem, fmin=fmin,
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
  }
  if (structof(x) != double) {
    x = double(x);
  }
  local gx, gnorm, isfree, iter, step;
  for (;;) {
    if (task == 1) {
      /* Evaluate function and gradient. */
      if (eval >= maxeval) {
        /* Too many function evaluations.  We restore the variables at the
           start of the line search which is a cheap way (no extra memory cost)
           to recover variables which should be nearly the best ones. */
        stop = 1n;
        msg = swrite(format="warning: too many function evaluations (%d)\n",
                     eval);
        opl_vmlmb_restore, x, fx, gx, ws;
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
    if (task == 2 && bounds != 0) {
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
    if (task > 2) {
      iter = opl_vmlmb_get_iterations(ws);
      gnorm = opl_vmlmb_get_gnorm(ws);
      if (task > 3) {
        /* Error or warning. */
        stop = 1n;
        msg = opl_vmlmb_get_msg(ws);
      } else if (gnorm <= gtest) {
        stop = 1n;
        msg = swrite(format="convergence (%s)\n", "gradient small enough");
      } else if (iter > maxiter) {
        stop = 1n;
        msg = swrite(format="warning: too many iterations (%d)\n", iter);
      }
    }
    if (verb && (stop || task > 2 && (iter % verb) == 0)) {
      if (eval == 1 && ! use_printer) {
        write, output, format="# Method %d (MEM=%d): %s\n#\n",
          0, mem, method_name;
        write, output, format="# %s\n# %s\n",
          "ITER  EVAL   CPU (ms)        FUNC               GNORM   STEPLEN",
          "---------------------------------------------------------------";
      }
      timer, elapsed;
      cpu = 1e3*(elapsed(1) - cpu_start);
      step = opl_vmlmb_get_step(ws);
      if (task <= 2) {
        iter = opl_vmlmb_get_iterations(ws);
        gnorm = opl_vmlmb_get_gnorm(ws);
      }
      if (use_printer) {
        printer, output, iter, eval, cpu, fx, gnorm, step, x, extra;
      } else {
        write, output, format=" %5d %5d %10.3f  %+-24.15e%-9.1e%-9.1e\n",
          iter, eval, cpu, fx, gnorm, step;
      }
      if (use_viewer) {
        viewer, x, extra;
      }
    }
    if (stop) {
      if (msg && (verb || (task != 3 && ! quiet))) {
        write, output, format="# %s\n", strtrim(msg, 2, blank=" \t\v\n\r");
      }
      return x;
    }

    /* Call optimizer. */
    task = opl_vmlmb_next(x, fx, gx, ws, isfree);
  }
}

/*---------------------------------------------------------------------------*/
/* PRIVATE ROUTINES */

func __opl_vmlmb_parse(ws)
{
  extern m, n, csave, isave, dsave;
  if (structof(ws) != pointer || numberof(ws) != 3) {
    error, "expecting VMLMB workspace";
  }
  eq_nocopy, csave, *ws(1);
  if (structof(csave) != char || numberof(csave) != 128) {
    error, "corrupted workspace (CSAVE)";
  }
  eq_nocopy, isave, *ws(2);
  if (structof(isave) != long || numberof(isave) != 11) {
    error, "corrupted workspace (ISAVE)";
  }
  m = isave(4);
  n = isave(5);
  eq_nocopy, dsave, *ws(3);
  if (structof(dsave) != double || numberof(dsave) != 28 + n + 2*m*(1 + n)) {
    error, "corrupted workspace (DSAVE)";
  }

}
errs2caller, __opl_vmlmb_parse;

extern __opl_vmlmb_setup;
/* PROTOTYPE
   int opl_vmlmb_setup(long n, long m,
                       double fatol, double frtol,
                       double sftol, double sgtol, double sxtol,
                       double delta, double epsilon,
                       char array csave, long array isave, double array dsave);
*/

extern __opl_vmlmb_next;
/* PROTOTYPE
   int opl_vmlmb_next(double array x, double array f, double array g,
                      pointer isfree, pointer h,
		      char array csave, long array isave, double array dsave);
*/

extern __opl_vmlmb_restart;
/* PROTOTYPE
   int opl_vmlmb_restart(char array csave, long array isave, double array dsave);
*/

extern __opl_vmlmb_restore;
/* PROTOTYPE
   int opl_vmlmb_restore(double array x, double array f, double array g,
                         char array csave, long array isave, double array dsave);
*/

extern __opl_vmlmb_set_fmin;
/* PROTOTYPE
   int opl_vmlmb_set_fmin(char array csave,
                          long array isave,
                          double array dsave,
                          double new_value,
                          double array old_value);
*/
extern __opl_vmlmb_get_fmin;
/* PROTOTYPE
   int opl_vmlmb_get_fmin(char array csave,
                         long array isave,
                         double array dsave,
                         double array value);
*/
