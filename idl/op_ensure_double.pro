function op_ensure_double, x
;+
; NAME:
;   op_ensure_double
;
;
; PURPOSE:
;   Make sure a variable is of type DOUBLE.
;
;
; CATEGORY:
;   OptimPack.
;
;
; CALLING SEQUENCE:
;   op_ensure_double(x)
;
;
; INPUTS:
;   X = IDL variable to check/fix.
;
;
; OPTIONAL INPUTS:
;   None.
;
;
; KEYWORD PARAMETERS:
;   None.
;
;
; OUTPUTS:
;   Returned value is TRUE if X cannot be converted to type DOUBLE,
;   FALSE otherwise.
;
;
; OPTIONAL OUTPUTS:
;   None.
;
;
; COMMON BLOCKS:
;   None.
;
;
; SIDE EFFECTS:
;   If X is not of type DOUBLE, it is converted to DOUBLE (if possible).
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
;   The following piece of code ensures that X is converted to double or
;   raise an error:
;     if op_ensure_double(x) then message, "bad data type for X"
;
;
; MODIFICATION HISTORY:
;   2003, Eric THIEBAUT.
;   $Id$
;   $Log$
;-
  type = SIZE(x, /TYPE)
  if type eq 5L then return, 0b
  if (type ge 1L and type le 4L) or (type ge 12L and type le 15L) then begin
    x = double(x)
    return, 0b
  end
  return, 1b
end
