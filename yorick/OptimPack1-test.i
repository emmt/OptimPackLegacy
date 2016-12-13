/*
 * OptimPack1-test.i --
 *
 * Various tests from MINPACK suite for the optimization routines in
 * OptimPack extension for Yorick.
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

require, "OptimPack1.i";

#if 0
#include "OptimPack1-test.i"
for (i=1;i<=18;++i) op_test_um, prob=i, method=0, verb=1, output="OptimPack1-test.out";

// Since the CPU time may change, you can compare the outputs
// with:
//   sed -e 's/^\( *[0-9]* *[0-9]*\) *[^ ]*/\1/'
//
#endif


local op_rosenbrock_nevals;
op_rosenbrock_nevals=0;

func op_test_rosenbrock(nil, start=, method=, ndirs=, frtol=)
/* DOCUMENT op_test_rosenbrock, ...;
     Test op_driver with Rosenbrock function.

   SEE ALSO: op_test_rosenbrock_func. */
{
  x = is_void(start) ? [0.0, 0.0] : start;
  return op_driver(op_test_rosenbrock_func,
                   (is_void(start) ? [0.0, 0.0] : start),
                   method=method, fmin=0.0, verb=1, ndirs=ndirs,
                   frtol=frtol);
}

func op_test_rosenbrock_func(x, &g)
{
  // Rosenbrock:
  //    f  = 100*(x2 - x1^2)^2 + (1 - x1)^2
  x1 = x(1);
  u = x(2) - x1*x1;
  v = 1.0 - x1;
  f = 100.0*u*u + v*v;
  g = [-400.0*u*x1 - 2.0*v, 200.0*u];
  ++op_rosenbrock_nevals;
  return f;
}

func op_test_quad(x, &g)
{
  u = (x - [1.0, 1.0, 1.0]);
  w = [3.0, 7.0, 2.0];
  g = 2.0*w*u;
  return sum(w*u*u);
}

func op_test_um(prob=, n=, method=, ndir=, verb=, factor=,
                maxiter=, maxeval=, frtol=, fatol=,
                sftol=, sgtol=, sxtol=,
                output=)
/* DOCUMENT op_test_um(...)
     Check various optimization methods for eighteen nonlinear unconstrained
     minimization problems.

   SEE ALSO: optim_vmlm, optim_cgmn, minpack1_umobj. */
{
  /* Output array. */
  local result;

  /* Output stream. */
  if (! is_void(output)) {
    if (structof(output) == string) {
      output = open(output, "a");
    } else if (typeof(output) != "text_stream") {
      error, "bad value for keyword OUTPUT";
    }
  }

  /* Check compatibility of arguments N and PROB. */
  if (prob == 1) {
    name = "Helical valley function.";
    if (is_void(n)) n = 3;
    else if (n != 3) error, "N must be 3 for problem #1";
  } else if (prob == 2) {
    name = "Biggs exp6 function.";
    if (is_void(n)) n = 6;
    else if (n != 6) error, "N must be 6 for problem #2";
  } else if (prob == 3) {
    name = "Gaussian function.";
    if (is_void(n)) n = 3;
    else if (n != 3) error, "N must be 3 for problem #3";
  } else if (prob == 4) {
    name = "Powell badly scaled function.";
    if (is_void(n)) n = 2;
    else if (n != 2) error, "N must be 2 for problem #4";
  } else if (prob == 5) {
    name = "Box 3-dimensional function.";
    if (is_void(n)) n = 3;
    else if (n != 3) error, "N must be 3 for problem #5";
  } else if (prob == 6) {
    name = "Variably dimensioned function.";
    if (is_void(n)) n = 10;
    else if (n < 1) error, "N must be >= 1 in problem #6";
  } else if (prob == 7) {
    name = "Watson function.";
    msg = "N may be 2 or greater but is usually 6 or 9 for problem #7";
    if (is_void(n)) {
      write, msg;
      n = 6;
    } else if (n < 2) error, msg;
  } else if (prob == 8) {
    name = "Penalty function I.";
    if (is_void(n)) n = 10;
    else if (n < 1) error, "N must be >= 1 in problem #8";
  } else if (prob == 9) {
    name = "Penalty function II.";
    if (is_void(n)) n = 10;
    else if (n < 1) error, "N must be >= 1 in problem #9";
  } else if (prob == 10) {
    name = "Brown badly scaled function.";
    if (is_void(n)) n = 2;
    else if (n != 2) error, "N must be 2 for problem #10";
  } else if (prob == 11) {
    name = "Brown and Dennis function.";
    if (is_void(n)) n = 4;
    else if (n != 4) error, "N must be 4 for problem #11";
  } else if (prob == 12) {
    name = "Gulf research and development function.";
    if (is_void(n)) n = 3;
    else if (n != 3) error, "N must be 3 for problem #12";
  } else if (prob == 13) {
    name = "Trigonometric function.";
    if (is_void(n)) n = 10;
    else if (n < 1) error, "N must be >= 1 in problem #13";
  } else if (prob == 14) {
    name = "Extended Rosenbrock function.";
    if (is_void(n)) n = 10;
    else if (n < 1 || n%2 != 0)
      error, "N must be a multiple of 2 in problem #14";
  } else if (prob == 15) {
    name = "Extended Powell function.";
    if (is_void(n)) n = 12;
    else if (n < 1 || n%4 != 0)
      error, "N must be a multiple of 4 in problem #15";
  } else if (prob == 16) {
    name = "Beale function.";
    if (is_void(n)) n = 2;
    else if (n != 2) error, "N must be 2 for problem #16";
  } else if (prob == 17) {
    name = "Wood function.";
    if (is_void(n)) n = 4;
    else if (n != 4) error, "N must be 4 for problem #17";
  } else if (prob == 18) {
    name = "Chebyquad function.";
    if (is_void(n)) n = 25;
    else if (n<1 || n>50) error, "N must be <=50 for problem #18";
  } else error, "PROB must be an integer between 1 and 18";

  /* Maximum number of iterations and function evaluations. */
  if (is_void(maxiter)) maxiter = 100*n;
  if (is_void(maxeval)) maxeval = 10*maxiter;

  /* Starting vector. */
  x = minpack1_umipt(n, prob, factor);
  dims = dimsof(x);

  /* Choose minimization method. */
  //if (is_void(frtol)) frtol = 1e-10;
  //if (is_void(fatol)) fatol = 1e-10;
  if (is_void(method)) method = 0;
  if (method == 0) {
    /* Variable metric. */
    if (is_void(ndir)) ndir = n;
    method_name = swrite(format="Limited Memory BFGS (VMLM with NDIR=%d)",
                         ndir);
    ws = op_vmlmb_setup(n, ndir, fmin=0.0,
                        fatol=fatol, frtol=frtol,
                        sftol=sftol, sgtol=sgtol, sxtol=sxtol);
  } else if (method < 0) {
    if (is_void(ndir)) ndir = n;
    method_name = swrite(format="Limited Memory BFGS (LBFGS with NDIR=%d)",
                         ndir);
    ws = op_lbfgs_setup(n, ndir);
  } else if (method >= 1 && method <= 15) {
    /* Conjugate gradient. */
    ndir = 2;
    method_name = swrite(format="Conjugate Gradient (%s)",
                         ["Fletcher-Reeves", "Polak-Ribiere",
                          "Polak-Ribiere with non-negative BETA"](method&3));
    ws = optim_cgmn_setup(method, fmin=fmin, fatol=fatol, frtol=frtol);
  } else {
    error, "bad METHOD";
  }
  step = 0.0;
  task = 1;
  eval = iter = 0;
  elapsed = array(double, 3);
  timer, elapsed;
  cpu_start = elapsed(1);
  for (;;) {
    if (task == 1) {
      /* Evaluate function. */
      f = minpack1_umobj(x, prob);
      g = minpack1_umgrd(x, prob);
      ++eval;
    }

    /* Check for convergence. */
    if (task != 1 || eval == 1) {
      too_many_eval = (maxeval >= 1 && eval > maxeval);
      too_many_iter = (maxiter >= 1 && iter > maxiter);
      timer, elapsed;
      cpu = elapsed(1) - cpu_start;
      gnorm = sqrt(sum(g*g));
      if (verb) {
        if (eval == 1) {
          write, output, format="#\n# Problem %d (N=%d): %s\n",
            prob, n, name;
          write, output, format="# Method %d (NDIR=%d): %s\n#\n",
            method, ndir, method_name;
          write, output, format="# %s\n# %s\n",
            "ITER  EVAL   CPU (ms)        FUNC               GNORM   STEPLEN",
            "---------------------------------------------------------------";
        }
        write, output, format=" %5d %5d %10.3f  %+-24.15e%-9.1e%-9.1e\n",
          iter, eval, 1e3*cpu, f, gnorm, step;
      }
      if (! am_subroutine()) grow, result, [[iter, eval, 1e3*cpu, f,
                                             gnorm, step]];
      if (task > 2) {
        msg = op_vmlmb_msg(ws);
        break;
      }
      if (maxiter >= 0 && iter > maxiter && task == 2) {
        msg = swrite(format="warning: too many iterations (%d)\n",
                     iter);
        break;
      }
      if (maxeval >= 0 && eval > maxeval) {
        x = x0;
        msg = swrite(format="warning: too many function evaluation (%d)\n",
                     eval);
        break;
      }
    }


    /* Call optimizer. */
    if (! method) {
      task = op_vmlmb_next(x, f, g, ws);
      iter = (*ws(2))(7);
      step = (*ws(3))(22);
    } else if (method < 0) {
      task = op_lbfgs_next(x, f, g, ws);
      if (task == 2 || task == 3) ++iter;
      step = -1.0;
    } else {
      optim_cgmn, x, f, g, task, ws;
      iter = optim_cgmn_iter(ws);
      step = optim_cgmn_step(ws);
    }
  }

  if (task == 3) status = 0;
  else if (ident == 4) status = 1;
  else status = 2;
  write, output, format=(verb?"# %s\n#\n":"*** %s\n"), msg;
  return result;
  //return [status, n, ndir, method, iter, eval, cpu, f];
  //return x;
}

/*---------------------------------------------------------------------------*/
/* NONLINEAR UNCONSTRAINED MINIMIZATION PROBLEMS
 *
 *   This suite of problems is taken from the MINPACK Project.
 *
 *   HISTORY:
 *     - Argonne National Laboratory. MINPACK Project. March 1980.
 *       Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More.
 *     - Conversion to Yorick. November 2001. Eric Thiebaut.
 */
_um_y = [9.0e-4,   4.4e-3,   1.75e-2,  5.4e-2,   1.295e-1,
         2.42e-1,  3.521e-1, 3.989e-1, 3.521e-1, 2.42e-1,
         1.295e-1, 5.4e-2,   1.75e-2,  4.4e-3,   9.0e-4];

func minpack1_umobj(x, prob)
/* DOCUMENT minpack1_umobj(x, prob)
     Returns  the objective functions  of eighteen  nonlinear unconstrained
     minimization problems.  X  is the parameter array: a  vector of length
     N, PROB is the problem number (a positive integer between 1 and 18).

     The  values  of  N  for  functions 1,2,3,4,5,10,11,12,16  and  17  are
     3,6,3,2,3,2,4,3,2 and 4, respectively.  For  function 7, n may be 2 or
     greater but is usually 6 or 9.  For functions 6,8,9,13,14,15 and 18, N
     may be variable,  however it must be even for  function 14, a multiple
     of 4 for function 15, and not greater than 50 for function 18.

   SEE ALSO: minpack1_umgrd, minpack1_umipt. */
{
  n = numberof(x);

  /* Function routine selector. */
  if (prob == 1) {
    /* Helical valley function. */
    tpi = 8.0*atan(1.0);
    if      (x(1) > 0.0) th = atan(x(2)/x(1))/tpi;
    else if (x(1) < 0.0) th = atan(x(2)/x(1))/tpi + 0.5;
    else                 th = (x(2) >= 0.0 ? 0.25 : -0.25);
    arg = x(1)*x(1) + x(2)*x(2);
    r = sqrt(arg);
    t = x(3) - 10.0*th;
    f = 100.0*(t*t + (r - 1.0)*(r - 1.0)) + x(3)*x(3);
    return f;
  } else if (prob == 2) {
    /* Biggs exp6 function. */
    f = 0.0;
    for (i=1 ; i<=13 ; i+=1) {
      d1 = double(i)/10.0;
      d2 = exp(-d1) - 5.0*exp(-10.0*d1) + 3.0*exp(-4.0*d1);
      s1 = exp(-d1*x(1));
      s2 = exp(-d1*x(2));
      s3 = exp(-d1*x(5));
      t = x(3)*s1 - x(4)*s2 + x(6)*s3 - d2;
      f += t*t;
    }
    return f;
  } else if (prob == 3) {
    /* Gaussian function. */
    f = 0.0;
    for (i=1 ; i<=15 ; i+=1) {
      d1 = 0.5*double(i-1);
      d2 = 3.5 - d1 - x(3);
      arg = -0.5*x(2)*d2*d2;
      r = exp(arg);
      t = x(1)*r - _um_y(i);
      f += t*t;
    }
    return f;
  } else if (prob == 4) {
    /* Powell badly scaled function. */
    t1 = 1e4*x(1)*x(2) - 1.0;
    s1 = exp(-x(1));
    s2 = exp(-x(2));
    t2 = s1 + s2 - 1.0001;
    f = t1*t1 + t2*t2;
    return f;
  } else if (prob == 5) {
    /* Box 3-dimensional function. */
    f = 0.0;
    for (i=1 ; i<=10 ; i+=1) {
      d1 = double(i);
      d2 = d1/10.0;
      s1 = exp(-d2*x(1));
      s2 = exp(-d2*x(2));
      s3 = exp(-d2) - exp(-d1);
      t = s1 - s2 - s3*x(3);
      f += t*t;
    }
    return f;
  } else if (prob == 6) {
    /* Variably dimensioned function. */
    t1 = 0.0;
    t2 = 0.0;
    for (j=1 ; j<=n ; j+=1) {
      t1 += double(j)*(x(j) - 1.0);
      t = x(j) - 1.0;
      t2 += t*t;
    }
    t = t1*t1;
    f = t2 + t*(1.0 + t);
    return f;
  } else if (prob == 7) {
    /* Watson function. */
    f = 0.0;
    for (i=1 ; i<=29 ; i+=1) {
      d1 = double(i)/29.0;
      s1 = 0.0;
      d2 = 1.0;
      for (j=2 ; j<=n ; j+=1) {
	s1 += double(j-1)*d2*x(j);
	d2 = d1*d2;
      }
      s2 = 0.0;
      d2 = 1.0;
      for (j=1 ; j<=n ; j+=1) {
	s2 += d2*x(j);
	d2 = d1*d2;
      }
      t = s1 - s2*s2 - 1.0;
      f += t*t;
    }
    t = x(1)*x(1);
    t1 = x(2) - t - 1.0;
    f += t + t1*t1;
    return f;
  } else if (prob == 8) {
    /* Penalty function I. */
    t1 = -0.25;
    t2 = 0.0;
    for (j=1 ; j<=n ; j+=1) {
      t1 += x(j)*x(j);
      t = x(j) - 1.0;
      t2 += t*t;
    }
    f = 1e-5*t2 + t1*t1;
    return f;
  } else if (prob == 9) {
    /* Penalty function II. */
    t1 = -1.0;
    t2 = 0.0;
    t3 = 0.0;
    d1 = exp(0.1);
    d2 = 1.0;
    s2 = 0.0; /* avoid compiler warning about `s2' used uninitialized */
    for (j=1 ; j<=n ; j+=1) {
      t1 += double(n-j+1)*x(j)*x(j);
      s1 = exp(x(j)/10.0);
      if (j > 1) {
	s3 = s1 + s2 - d2*(d1 + 1.0);
	t2 += s3*s3;
	t = (s1 - 1.0/d1);
	t3 += t*t;
      }
      s2 = s1;
      d2 = d1*d2;
    }
    t = (x(1) - 0.2);
    f = 1e-5*(t2 + t3) + t1*t1 + t*t;
    return f;
  } else if (prob == 10) {
    /* Brown badly scaled function. */
    t1 = x(1) - 1e6;
    t2 = x(2) - 2e-6;
    t3 = x(1)*x(2) - 2.0;
    f = t1*t1 + t2*t2 + t3*t3;
    return f;
  } else if (prob == 11) {
    /* Brown and Dennis function. */
    f = 0.0;
    for (i=1 ; i<=20 ; i+=1) {
      d1 = double(i)/5.0;
      d2 = sin(d1);
      t1 = x(1) + d1*x(2) - exp(d1);
      t2 = x(3) + d2*x(4) - cos(d1);
      t = t1*t1 + t2*t2;
      f += t*t;
    }
    return f;
  } else if (prob == 12) {
    /* Gulf research and development function. */
    f = 0.0;
    d1 = 2.0/3.0;
    for (i=1 ; i<=99 ; i+=1) {
      arg = double(i)/100.0;
      r = (-50.0*log(arg))^d1 + 25.0 - x(2);
      t1 = (abs(r)^x(3))/x(1);
      t2 = exp(-t1);
      t = t2 - arg;
      f += t*t;
    }
    return f;
  } else if (prob == 13) {
    /* Trigonometric function. */
    s1 = 0.0;
    for (j=1 ; j<=n ; j+=1) {
      s1 += cos(x(j));
    }
    f = 0.0;
    for (j=1 ; j<=n ; j+=1) {
      t = double(n+j) - sin(x(j)) - s1 - double(j)*cos(x(j));
      f += t*t;
    }
    return f;
  } else if (prob == 14) {
    /* Extended Rosenbrock function. */
    f = 0.0;
    for (j=1 ; j<=n ; j+=2) {
      t1 = 1.0 - x(j);
      t2 = 10.0*(x(j+1) - x(j)*x(j));
      f += t1*t1 + t2*t2;
    }
    return f;
  } else if (prob == 15) {
    /* Extended Powell function. */
    f = 0.0;
    for (j=1 ; j<=n ; j+=4) {
      t = x(j) + 10.0*x(j+1);
      t1 = x(j+2) - x(j+3);
      s1 = 5.0*t1;
      t2 = x(j+1) - 2.0*x(j+2);
      s2 = t2*t2*t2;
      t3 = x(j) - x(j+3);
      s3 = 10.0*t3*t3*t3;
      f += t*t + s1*t1 + s2*t2 + s3*t3;
    }
    return f;
  } else if (prob == 16) {
    /* Beale function. */
    s1 = 1.0 - x(2);
    t1 = 1.5 - x(1)*s1;
    s2 = 1.0 - x(2)*x(2);
    t2 = 2.25 - x(1)*s2;
    s3 = 1.0 - x(2)*x(2)*x(2);
    t3 = 2.625 - x(1)*s3;
    f = t1*t1 + t2*t2 + t3*t3;
    return f;
  } else if (prob == 17) {
    /* Wood function. */
    s1 = x(2) - x(1)*x(1);
    s2 = 1.0 - x(1);
    s3 = x(2) - 1.0;
    t1 = x(4) - x(3)*x(3);
    t2 = 1.0 - x(3);
    t3 = x(4) - 1.0;
    f = 100.0*s1*s1 + s2*s2 + 90.0*t1*t1 + t2*t2 + \
      10.0*(s3 + t3)*(s3 + t3) + (s3 - t3)*(s3 - t3)/10.0;
    return f;
  } else if (prob == 18) {
    /* Chebyquad function. */
    fvec = array(0.0, n);
    for (j=1 ; j<=n ; j+=1) {
      t1 = 1.0;
      t2 = 2.0*x(j) - 1.0;
      t = 2.0*t2;
      for (i=1 ; i<=n ; i+=1) {
	fvec(i) += t2;
	th = t*t2 - t1;
	t1 = t2;
	t2 = th;
      }
    }
    f = 0.0;
    d1 = 1.0/double(n);
    iev = -1;
    for (i=1 ; i<=n ; i+=1) {
      t = d1*fvec(i);
      if (iev > 0) t += 1.0/(double(i)*double(i) - 1.0);
      f += t*t;
      iev = -iev;
    }
    return f;
  }
  return 0.0;
}

func minpack1_umgrd(x, prob)
/* DOCUMENT minpack1_umgrd(x, prob)
     Returns  the  gradient  vectors  of eighteen  nonlinear  unconstrained
     minimization problems. The problem  dimensions are as described in the
     prologue comments of minpack1_umobj.

   SEE ALSO: minpack1_umobj, minpack1_umipt. */
{
  n = numberof(x);
  g = array(double, n);

  /* Gradient routine selector. */
  if (prob == 1) {
    /* Helical valley function. */
    tpi = 8.0*atan(1.0);
    if      (x(1) > 0.0) th = atan(x(2)/x(1))/tpi;
    else if (x(1) < 0.0) th = atan(x(2)/x(1))/tpi + 0.5;
    else                 th = (x(2) >= 0.0 ? 0.25 : -0.25);
    arg = x(1)*x(1) + x(2)*x(2);
    r = sqrt(arg);
    t = x(3) - 10.0*th;
    s1 = 10.0*t/(tpi*arg);
    g(1) = 200.0*(x(1) - x(1)/r + x(2)*s1);
    g(2) = 200.0*(x(2) - x(2)/r - x(1)*s1);
    g(3) = 2.0*(100.0*t + x(3));
  } else if (prob == 2) {
    /* Biggs exp6 function. */
    for (j=1 ; j<=6 ; ++j) g(j) = 0.0;
    for (i=1 ; i<=13 ; ++i) {
      d1 = double(i)/10.0;
      d2 = exp(-d1) - 5.0*exp(-10.0*d1) + 3.0*exp(-4.0*d1);
      s1 = exp(-d1*x(1));
      s2 = exp(-d1*x(2));
      s3 = exp(-d1*x(5));
      t = x(3)*s1 - x(4)*s2 + x(6)*s3 - d2;
      th = d1*t;
      g(1) = g(1) - s1*th;
      g(2) = g(2) + s2*th;
      g(3) = g(3) + s1*t;
      g(4) = g(4) - s2*t;
      g(5) = g(5) - s3*th;
      g(6) = g(6) + s3*t;
    }
    g(1) = 2.0*x(3)*g(1);
    g(2) = 2.0*x(4)*g(2);
    g(3) = 2.0*g(3);
    g(4) = 2.0*g(4);
    g(5) = 2.0*x(6)*g(5);
    g(6) = 2.0*g(6);
  } else if (prob == 3) {
    /* Gaussian function. */
    g(1) = 0.0;
    g(2) = 0.0;
    g(3) = 0.0;
    for (i=1 ; i<=15 ; ++i) {
      d1 = 0.5*double(i-1);
      d2 = 3.5 - d1 - x(3);
      arg = -0.5*x(2)*d2*d2;
      r = exp(arg);
      t = x(1)*r - _um_y(i);
      s1 = r*t;
      s2 = d2*s1;
      g(1) = g(1) + s1;
      g(2) = g(2) - d2*s2;
      g(3) = g(3) + s2;
    }
    g(1) = 2.0*g(1);
    g(2) = x(1)*g(2);
    g(3) = 2.0*x(1)*x(2)*g(3);
  } else if (prob == 4) {
    /* Powell badly scaled function. */
    t1 = 1e4*x(1)*x(2) - 1.0;
    s1 = exp(-x(1));
    s2 = exp(-x(2));
    t2 = s1 + s2 - 1.0001;
    g(1) = 2.0*(1e4*x(2)*t1 - s1*t2);
    g(2) = 2.0*(1e4*x(1)*t1 - s2*t2);
  } else if (prob == 5) {
    /* Box 3-dimensional function. */
    g(1) = 0.0;
    g(2) = 0.0;
    g(3) = 0.0;
    for (i=1 ; i<=10 ; ++i) {
      d1 = double(i);
      d2 = d1/10.0;
      s1 = exp(-d2*x(1));
      s2 = exp(-d2*x(2));
      s3 = exp(-d2) - exp(-d1);
      t = s1 - s2 - s3*x(3);
      th = d2*t;
      g(1) = g(1) - s1*th;
      g(2) = g(2) + s2*th;
      g(3) = g(3) - s3*t;
    }
    g(1) = 2.0*g(1);
    g(2) = 2.0*g(2);
    g(3) = 2.0*g(3);
  } else if (prob == 6) {
    /* Variably dimensioned function. */
    t1 = 0.0;
    for (j=1 ; j<=n ; ++j) {
      t1 += double(j)*(x(j) - 1.0);
    }
    t = t1*(1.0 + 2.0*t1*t1);
    for (j=1 ; j<=n ; ++j) {
      g(j) = 2.0*(x(j) - 1.0 + double(j)*t);
    }
  } else if (prob == 7) {
    /* Watson function. */
    for (j=1 ; j<=n ; ++j) {
      g(j) = 0.0;
    }
    for (i=1 ; i<=29 ; ++i) {
      d1 = double(i)/29.0;
      s1 = 0.0;
      d2 = 1.0;
      for (j=2 ; j<=n ; ++j) {
	s1 += double(j-1)*d2*x(j);
	d2 = d1*d2;
      }
      s2 = 0.0;
      d2 = 1.0;
      for (j=1 ; j<=n ; ++j) {
	s2 += d2*x(j);
	d2 = d1*d2;
      }
      t = s1 - s2*s2 - 1.0;
      s3 = 2.0*d1*s2;
      d2 = 2.0/d1;
      for (j=1 ; j<=n ; ++j) {
	g(j) = g(j) + d2*(double(j-1) - s3)*t;
	d2 = d1*d2;
      }
    }
    t1 = x(2) - x(1)*x(1) - 1.0;
    g(1) = g(1) + x(1)*(2.0 - 4.0*t1);
    g(2) = g(2) + 2.0*t1;
  } else if (prob == 8) {
    /* Penalty function I. */
    t1 = -0.25;
    for (j=1 ; j<=n ; ++j) {
      t1 += x(j)*x(j);
    }
    d1 = 2.0*1e-5;
    th = 4.0*t1;
    for (j=1 ; j<=n ; ++j) {
      g(j) = d1*(x(j) - 1.0) + x(j)*th;
    }
  } else if (prob == 9) {
    /* Penalty function II. */
    s2 = 0.0; /* avoid compiler warning about `s2' used uninitialized */
    t1 = -1.0;
    for (j=1 ; j<=n ; ++j) {
      t1 += double(n-j+1)*x(j)*x(j);
    }
    d1 = exp(0.1);
    d2 = 1.0;
    th = 4.0*t1;
    for (j=1 ; j<=n ; ++j) {
      g(j) = double(n-j+1)*x(j)*th;
      s1 = exp(x(j)/10.0);
      if (j > 1) {
	s3 = s1 + s2 - d2*(d1 + 1.0);
	g(j) = g(j) + 1e-5*s1*(s3 + s1 - 1.0/d1)/5.0;
	g(j-1) = g(j-1) + 1e-5*s2*s3/5.0;
      }
      s2 = s1;
      d2 = d1*d2;
    }
    g(1) = g(1) + 2.0*(x(1) - 0.2);
  } else if (prob == 10) {
    /* Brown badly scaled function. */
    t1 = x(1) - 1e6;
    t2 = x(2) - 2e-6;
    t3 = x(1)*x(2) - 2.0;
    g(1) = 2.0*(t1 + x(2)*t3);
    g(2) = 2.0*(t2 + x(1)*t3);
  } else if (prob == 11) {
    /* Brown and Dennis function. */
    g(1) = 0.0;
    g(2) = 0.0;
    g(3) = 0.0;
    g(4) = 0.0;
    for (i=1 ; i<=20 ; ++i) {
      d1 = double(i)/5.0;
      d2 = sin(d1);
      t1 = x(1) + d1*x(2) - exp(d1);
      t2 = x(3) + d2*x(4) - cos(d1);
      t = t1*t1 + t2*t2;
      s1 = t1*t;
      s2 = t2*t;
      g(1) = g(1) + s1;
      g(2) = g(2) + d1*s1;
      g(3) = g(3) + s2;
      g(4) = g(4) + d2*s2;
    }
    g(1) = 4.0*g(1);
    g(2) = 4.0*g(2);
    g(3) = 4.0*g(3);
    g(4) = 4.0*g(4);
  } else if (prob == 12) {
    /* Gulf research and development function. */
    g(1) = 0.0;
    g(2) = 0.0;
    g(3) = 0.0;
    d1 = 2.0/3.0;
    for (i=1 ; i<=99 ; ++i) {
      arg = double(i)/100.0;
      r = (-50.0*log(arg))^d1 + 25.0 - x(2);
      t1 = (abs(r)^x(3))/x(1);
      t2 = exp(-t1);
      t = t2 - arg;
      s1 = t1*t2*t;
      g(1) = g(1) + s1;
      g(2) = g(2) + s1/r;
      g(3) = g(3) - s1*log(abs(r));
    }
    g(1) = 2.0*g(1)/x(1);
    g(2) = 2.0*x(3)*g(2);
    g(3) = 2.0*g(3);
  } else if (prob == 13) {
    /* Trigonometric function. */
    s1 = 0.0;
    for (j=1 ; j<=n ; ++j) {
      g(j) = cos(x(j));
      s1 += g(j);
    }
    s2 = 0.0;
    for (j=1 ; j<=n ; ++j) {
      th = sin(x(j));
      t = double(n+j) - th - s1 - double(j)*g(j);
      s2 += t;
      g(j) = (double(j)*th - g(j))*t;
    }
    for (j=1 ; j<=n ; ++j) {
      g(j) = 2.0*(g(j) + sin(x(j))*s2);
    }
  } else if (prob == 14) {
    /* Extended Rosenbrock function. */
    for (j=1 ; j<=n ; j+=2) {
      t1 = 1.0 - x(j);
      g(j+1) = 200.0*(x(j+1) - x(j)*x(j));
      g(j) = -2.0*(x(j)*g(j+1) + t1);
    }
  } else if (prob == 15) {
    /* Extended Powell function. */
    for (j=1 ; j<=n ; j+=4) {
      t = x(j) + 10.0*x(j+1);
      t1 = x(j+2) - x(j+3);
      s1 = 5.0*t1;
      t2 = x(j+1) - 2.0*x(j+2);
      s2 = 4.0*t2*t2*t2;
      t3 = x(j) - x(j+3);
      s3 = 20.0*t3*t3*t3;
      g(j) = 2.0*(t + s3);
      g(j+1) = 20.0*t + s2;
      g(j+2) = 2.0*(s1 - s2);
      g(j+3) = -2.0*(s1 + s3);
    }
  } else if (prob == 16) {
    /* Beale function. */
    s1 = 1.0 - x(2);
    t1 = 1.5 - x(1)*s1;
    s2 = 1.0 - x(2)*x(2);
    t2 = 2.25 - x(1)*s2;
    s3 = 1.0 - x(2)*x(2)*x(2);
    t3 = 2.625 - x(1)*s3;
    g(1) = -2.0*(s1*t1 + s2*t2 + s3*t3);
    g(2) = 2.0*x(1)*(t1 + x(2)*(2.0*t2 + 3.0*x(2)*t3));
  } else if (prob == 17) {
    /* Wood function. */
    s1 = x(2) - x(1)*x(1);
    s2 = 1.0 - x(1);
    s3 = x(2) - 1.0;
    t1 = x(4) - x(3)*x(3);
    t2 = 1.0 - x(3);
    t3 = x(4) - 1.0;
    g(1) = -2.0*(200.0*x(1)*s1 + s2);
    g(2) = 200.0*s1 + 20.2*s3 + 19.8*t3;
    g(3) = -2.0*(180.0*x(3)*t1 + t2);
    g(4) = 180.0*t1 + 20.2*t3 + 19.8*s3;
  } else if (prob == 18) {
    /* Chebyquad function. */
    fvec = array(0.0, n);
    for (j=1 ; j<=n ; ++j) {
      t1 = 1.0;
      t2 = 2.0*x(j) - 1.0;
      t = 2.0*t2;
      for (i=1 ; i<=n ; ++i) {
	fvec(i) += t2;
	th = t*t2 - t1;
	t1 = t2;
	t2 = th;
      }
    }
    d1 = 1.0/double(n);
    iev = -1;
    for (i=1 ; i<=n ; ++i) {
      fvec(i) *= d1;
      if (iev > 0) fvec(i) += 1.0/(double(i)*double(i) - 1.0);
      iev = -iev;
    }
    for (j=1 ; j<=n ; ++j) {
      g(j) = 0.0;
      t1 = 1.0;
      t2 = 2.0*x(j) - 1.0;
      t = 2.0*t2;
      s1 = 0.0;
      s2 = 2.0;
      for (i=1 ; i<=n ; ++i) {
	g(j) = g(j) + fvec(i)*s2;
	th = 4.0*t2 + t*s2 - s1;
	s1 = s2;
	s2 = th;
	th = t*t2 - t1;
	t1 = t2;
	t2 = th;
      }
    }
    d2 = 2.0*d1;
    for (j=1 ; j<=n ; ++j) g(j) *= d2;
  }
  return g;
}

func minpack1_umipt(n, prob, factor)
/* DOCUMENT minpack1_umipt(n, prob, factor)
     Returns  the standard  starting points  for the  functions  defined by
     subroutine  minpack1_umobj.  The  function  returns a  vector  X of  N
     elements, X is a multiple  (times FACTOR, default 1.0) of the standard
     starting point.  For the  seventh function the standard starting point
     is 0.0,  so in this  case, if FACTOR  is not unity, then  the function
     returns  X  filled with  FACTOR.   PROB has  the  same  meaning as  in
     minpack1_umobj.

   SEE ALSO: minpack1_umobj, minpack1_umipt. */
{
  if (is_void(factor)) factor = 1.0;
  x = array(double, n);

  /* Selection of initial point. */
  if (prob == 1) {
    /* Helical valley function. */
    x(1) = -1.0;
    x(2) = 0.0;
    x(3) = 0.0;
  } else if (prob == 2) {
    /* Biggs exp6 function. */
    x(1) = 1.0;
    x(2) = 2.0;
    x(3) = 1.0;
    x(4) = 1.0;
    x(5) = 1.0;
    x(6) = 1.0;
  } else if (prob == 3) {
    /* Gaussian function. */
    x(1) = 0.4;
    x(2) = 1.0;
    x(3) = 0.0;
  } else if (prob == 4) {
    /* Powell badly scaled function. */
    x(1) = 0.0;
    x(2) = 1.0;
  } else if (prob == 5) {
    /* Box 3-dimensional function. */
    x(1) = 0.0;
    x(2) = 10.0;
    x(3) = 20.0;
  } else if (prob == 6) {
    /* Variably dimensioned function. */
    h = 1.0/double(n);
    for (j=1 ; j<=n ; ++j) x(j) = 1.0 - double(j)*h;
  } else if (prob == 7) {
    /* Watson function. */
    for (j=1 ; j<=n ; ++j) x(j) = 0.0;
  } else if (prob == 8) {
    /* Penalty function I. */
    for (j=1 ; j<=n ; ++j) x(j) = double(j);
  } else if (prob == 9) {
    /* Penalty function II. */
    for (j=1 ; j<=n ; ++j) x(j) = 0.5;
  } else if (prob == 10) {
    /* Brown badly scaled function. */
    x(1) = 1.0;
    x(2) = 1.0;
  } else if (prob == 11) {
    /* Brown and Dennis function. */
    x(1) = 25.0;
    x(2) = 5.0;
    x(3) = -5.0;
    x(4) = -1.0;
  } else if (prob == 12) {
    /* Gulf research and development function. */
    x(1) = 5.0;
    x(2) = 2.5;
    x(3) = 0.15;
  } else if (prob == 13) {
    /* Trigonometric function. */
    h = 1.0/double(n);
    for (j=1 ; j<=n ; ++j) x(j) = h;
  } else if (prob == 14) {
    /* Extended Rosenbrock function. */
    for (j=1 ; j<=n ; j+=2) {
      x(j) = -1.2;
      x(j+1) = 1.0;
    }
  } else if (prob == 15) {
    /* Extended Powell singular function. */
    for (j=1 ; j<=n ; j+=4) {
      x(j) = 3.0;
      x(j+1) = -1.0;
      x(j+2) = 0.0;
      x(j+3) = 1.0;
    }
  } else if (prob == 16) {
    /* Beale function. */
    x(1) = 1.0;
    x(2) = 1.0;
  } else if (prob == 17) {
    /* Wood function. */
    x(1) = -3.0;
    x(2) = -1.0;
    x(3) = -3.0;
    x(4) = -1.0;
  } else if (prob == 18) {
    /* Chebyquad function. */
    h = 1.0/double(n+1);
    for (j=1 ; j<=n ; ++j) x(j) = double(j)*h;
  }

  /* Compute multiple of initial point. */
  if (factor != 1.0) {
    if (prob == 7) {
      for (j=1 ; j<=n ; ++j) x(j) = factor;
    } else {
      for (j=1 ; j<=n ; ++j) x(j) *= factor;
    }
  }
  return x;
}
