function op_vmlmb, x, f, g, ws, active=active, h=h
;+
; NAME:
;   op_vmlmb
;
;
; PURPOSE:
;   Computes a local minimizer of a function of several variables by a
;   limited memory variable metric (BFGS) method; optionally, the
;   parameters may be bounded.
;
;
; CATEGORY:
;   OptimPack.
;
;
; CALLING SEQUENCE:
;   job = op_vmlmb(x, f, g, ws[, active[, h]])
;
;
; INPUTS:
;   X is a variable that contains the parameters (an array of length N,
;       see WS below).  On entry, X is an approximation to the solution.
;       On exit with JOB=3, X is the current approximation.
;
;   F is a scalar variable.  On entry, F is the value of the function at X.
;       On final exit, F is the function value at X.
;
;   G is a variable that contains the gradient (an array of length N).  On
;       entry, G is the value of the gradient at X.  On final exit, G is
;       the value of the gradient at the solution X.
;
;   WS is a workspace array as returned by op_vmlmb_setup.  WS must not
;       be changed between calls to op_vmlmb and the value of N when WS
;       was created must match be the number of elements of X (and G).
;
;
; OPTIONAL INPUTS:
;   None.
;
;
; KEYWORD PARAMETERS:
;   ACTIVE is an optional logical array with length N provided by the
;       caller if the values in X has bounds.  If the parameters have no
;       bounds, ACTIVE should be undefined or not specified (unconstrained
;       minimization).  Otherwise, elements set to zero (false) in ACTIVE
;       indicate that the corresponding values in X has reached a bound and
;       should not be changed during the next step because the gradient has
;       the wrong sign (i.e.  the steepest descent direction would violate
;       the bound constraints):
;           ACTIVE[i] = 0 if X[i] has a lower bound XLO[i]
;                           and X[i]=XLO[i] and G[i]>=0 
;                       0 if X[i] value has an upper bound XHI[i]
;                           and X[i]=XHI[i] and G[i]<=0
;                       1 (or any non-zero value) otherwise
;
;       ACTIVE needs only to be computed (and specified) the first time
;       op_vmlmb is called and when JOB=2 (i.e.  after a successful step).
;       ACTIVE may also be specified when TASK=3 (i.e.  after convergence
;       if caller wish to continue with minimization).  If X has (some)
;       bounds, the caller is responsible for applying the bounds to X
;       before evaluating the function value F and the gradient G
;       (i.e. when JOB=1), e.g.:
;           if (X[i] < XLO[i]) X[i] = XLO[i];
;           if (X[i] > XHI[i]) X[i] = XHI[i];
;
;       If H is not specified or undefined or if H[i] > 0 for all i such
;       that ACTIVE[i] is non-zero, then ACTIVE is left unchanged.
;
;   H is an optional array with length N provided by the caller and such
;       that diag(H) is an approximation of the inverse of the Hessian
;       matrix.  If H is not specified or is undefined, then the inverse of
;       the Hessian is approximated by a simple rescaling using Shanno &
;       Phua formula.  Otherwise, if ACTIVE is undefined, all elements of H
;       must be strictly greater than zero; else ACTIVE[i] is set to zero
;       (false) if H[i] <= 0 (this is the only case where ACTIVE is
;       modified).  As for ACTIVE, H needs only to be specifed the first
;       time op_vmlmb is called and when JOB=2.
;
;
; OUTPUTS:
;   The returned value JOB is a scalar integer.  It can have one of the
;   following values:
;     JOB=1 - caller must evaluate the function and gradient at X and call
;             op_vmlm.
;     JOB=2 - a new iterate has been computed.  The approximation X,
;             function F, and gradient G are available for examination.
;     JOB=3 - the  search  is  successful.   The  solution,
;             function value and gradient are available in X, F and G.
;     JOB=4 - VMLMB is not able to satisfy the convergence conditions.  The
;             exit value of X contains the best approximation found so far.
;             Warning message can be obtained by op_vmlmb_msg(ws).
;     JOB=5 - there is an error in the input arguments.  Error message can
;             be obtained by op_vmlmb_msg(ws).
;
;
; OPTIONAL OUTPUTS:
;   None.
;
;
; COMMON BLOCKS:
;   OP_COMMON - used to store the name of the shared OptimPack library.
;
;
; SIDE EFFECTS:
;   X, F, G, ACTIVE and H may be converted to appropriate data types and
;   X, F and G are updated as explained above.
;
;
; RESTRICTIONS:
;   None.
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;
;   A typical invocation of VMLMB for unconstrained minimization has the
;   following outline:
;
;    ;; Choose a starting solution:
;    x = ...
;  
;    ;; Allocate and setup workspaces:
;    ws = op_vmlmb_setup(n_elements(x), ...)
;    job = 1
;
;    ;; Optimization loop:
;    while 1 do begin
;      if (job eq 1) begin
;        ;; Evaluate the function and the gradient at X
;        f = ...
;        g = ...
;      end else if (job eq 2) begin
;        ;; New successful step: the approximation X, function F, and
;        ;; gradient G, are available for inspection.
;        ...
;      end else begin
;        ;; Convergence (JOB=3), or warning (JOB=4), or error
;        print, op_vmlmb_msg(ws)
;        break
;      end
;      ;; Computes next step:
;      job = op_vmlmb(x, f, g, ws)
;    end
;
; A typical invocation of VMLMB for bound-constrained minimization has the
; following outline (XMIN and XMAX are the lower and upper bounds for X,
; the bounds may be arrays or scalar and it is easy to modify the example
; below if there is only an upper or a lower bound or if only some
; parameters are bounded):
;
;    ;; Choose a starting solution:
;    x = ...
;  
;    ;; Allocate and setup workspaces:
;    ws = op_vmlmb_setup(n_elements(x), ...)
;    job = 1
;
;    ;; Optimization loop:
;    eval = 0L     ; number of evaluations
;    while 1 do begin
;      if (job eq 1) begin
;        ;; Aply bound constraints
;        x = (x > xmin) < xmax
;
;        ;; Evaluate the function and the gradient at X
;        f = ...
;        g = ...
;        eval = eval + 1L
;      end else if (job eq 2) begin
;        ;; New successful step: the approximation X, function F, and
;        ;; gradient G, are available for inspection.
;        ...
;      end else begin
;        ;; Convergence (JOB=3), or warning (JOB=4), or error
;        print, op_vmlmb_msg(ws)
;        break
;      end
;      ;; Computes next step:
;      if (eval eq 1L or job ne 1) begin
;        ;; Computes set of active parameters:
;        active = (x gt xmin or g lt 0d0) and (x lt xmax or g gt 0d0)
;      end
;      job = op_vmlmb(x, f, g, ws, active=active)
;    end
;
;
;
; MODIFICATION HISTORY:
;   2003, Eric THIEBAUT.
;   $Id: op_vmlmb.pro,v 1.1 2007/07/11 05:50:35 eric Exp $
;   $Log: op_vmlmb.pro,v $
;   Revision 1.1  2007/07/11 05:50:35  eric
;   Initial revision
;
;-
  common op_common, libname
  on_error, 2
  if op_ensure_double(x) then message, "bad data type for X"
  if op_ensure_double(f) then message, "bad data type for F"
  if op_ensure_double(g) then message, "bad data type for G"
  type = SIZE(active, /TYPE)
  if arg_present(h) then begin
    if op_ensure_double(h) then message, "bad data type for H"
    if type ne 1L then begin
      ;; fix ACTIVE array
      if type eq 0L then active = h gt 0.0d0 else active = active ne 0
    end
    job = call_external(libname, 'op_idl_vmlmb_next', $
                        size(x), x, size(f), f, size(g), g, size(ws), ws, $
                        size(active), active, size(h), h)
  end else if type ne 0L then begin
    if type ne 1L then active = active ne 0
    job = call_external(libname, 'op_idl_vmlmb_next', $
                        size(x), x, size(f), f, size(g), g, size(ws), ws, $
                        size(active), active)
  end else begin
    job = call_external(libname, 'op_idl_vmlmb_next', $
                        size(x), x, size(f), f, size(g), g, size(ws), ws)
  end
  if job lt 0L then message, op_last_error()
  return, job
end
