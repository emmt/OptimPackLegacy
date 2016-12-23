function op_vmlmb_setup, n, m, fmin=fmin, frtol=frtol, fatol=fatol, $
                         sftol=sftol, sgtol=sgtol, sxtol=sxtol, $
                         epsilon=epsilon, delta=delta
;+
; NAME:
;   op_vmlmb_setup
;
;
; PURPOSE:
;   Set up workspace for op_vmlmb routine.
;
;
; CATEGORY:
;   OptimPack.
;
;
; CALLING SEQUENCE:
;   op_vmlmb_setup(n [, m])
;
;
; INPUTS:
;   N = Number of parameters.
;
;
; OPTIONAL INPUTS:
;   M = Number of memorized correction pairs used to compute the limited
;       memory variable metric (BFGS) approximation of the inverse of the
;       Hessian.  For large problems, M=3 to 5 gives good results.  For
;       small problems, M should be less or equal N.  The larger is M (and
;       N) the more computer memory will be needed to store the
;       workspaces. Default value is M=5.
;
;
; KEYWORD PARAMETERS:
;   FRTOL = Relative error desired in the function.  Convergence occurs if
;       the estimate of the relative error between F(X) and F(XSOL), where
;       XSOL is a local minimizer, is less or equal FRTOL.  FRTOL must have
;       a non-negative floating point value.  Default value is FATOL=1D-10.
;
;   FATOL = Absolute error desired in the function.  Convergence occurs if
;       the estimate of the absolute error between F(X) and F(XSOL), where
;       XSOL is a local minimizer, is less or equal FATOL.  FATOL must have
;       a non-negative floating point value.  Default value is FATOL=1D-13.
;
;   FMIN = Lower bound for the function.  VMLMB exits with a warning if
;       F < FMIN.  Default value is FMIN=0.0D0.
;
;   SFTOL, SGTOL, and SXTOL are tolerances for the line search subroutine
;       (see op_csrch).  Default values: SFTOL=0.001, SGTOL=0.9, and
;       SXTOL=0.1 (other values may be more suitable for highly
;       non-quadratic penalty function).
;
;
; OUTPUTS:
;   Returned value is a worspace array of type DOUBLE.
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
;   None.
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
;
;
; MODIFICATION HISTORY:
;   2003, Eric THIEBAUT.
;   Initial revision
;
;   2016, Eric THIEBAUT.
;   Update to OptimPackLegacy 1.4
;-
  common op_common, libname
  on_error, 2

  ;; Provides default values for parameters:
  if not n_elements(m)        then m = (5L < n)
  if not n_elements(fmin)     then fmin = 0.0D0
  if not n_elements(frtol)    then frtol = 1.0D-10
  if not n_elements(fatol)    then fatol = 1.0D-13
  if not n_elements(sftol)    then sftol = 1.0D-3
  if not n_elements(sgtol)    then sgtol = 9.0D-1
  if not n_elements(sxtol)    then sxtol = 1.0D-1
  if not n_elements(epsilon)  then epsilon = 0.0D0
  if not n_elements(delta)    then delta = 1.0D-2

  ;; Create workspace array (round the number of needed bytes up to the
  ;; size of a double)
  nbytes = call_external(libname, 'op_idl_vmlmb_size', $
                         size(n), n, size(m), m)
  ndoubles = ((nbytes + double_size - 1L)/double_size)
  ws = dblarr(ndoubles, /nozero)

  ;; Instanciate workspace array.
  task = call_external(libname, 'op_idl_vmlmb_init', $
                       size(ws),      ws,      $
                       size(n),       n,       $
                       size(m),       m,       $
                       size(fatol),   fatol,   $
                       size(frtol),   frtol,   $
                       size(delta),   delta,   $
                       size(epsilon), epsilon, $
                       size(sftol),   sftol,   $
                       size(sgtol),   sgtol,   $
                       size(sxtol),   sxtol)
  if task ne 1L then begin
    if task lt 0L then msg = op_last_error() else msg = op_vmlmb_msg(ws)
    message, msg
  endif
  return, ws
end
