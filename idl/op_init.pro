pro op_init, name
;+
; NAME:
;   op_init
;
;
; PURPOSE:
;   Initialize OptimPack routines.
;
;
; CATEGORY:
;   OptimPack.
;
;
; CALLING SEQUENCE:
;   op_init[, name]
;
;
; INPUTS:
;   NAME - Name (without extension) of the OptimPack-IDL library. Default
;          is "OptimPack_IDL".
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
;   None.
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
;   Data in common block OP_COMMON is updated.
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
;   op_init, "OptimPack_IDL"
;
;
; MODIFICATION HISTORY:
;   2003, Eric THIEBAUT.
;   $Id$
;   $Log$
;-
  common op_common, libname
  on_error, 2
  if n_elements(name) eq 0L then begin
    basename = 'OptimPack_IDL'
  end else begin
    basename = name
  end
  prefix = ''
  case !version.os_family of
    'unix': begin
      case !version.os of
        'hp-ux': suffix = '.sl'
        'AIX': begin
          ;; AIX won't find a shared lib in the current dir
          ;; unless the name is preceded with a ./
          suffix = '.a'
          if strmid(basename, 0, 1) ne '/' then prefix = './'
        end
        'sunos': begin
          suffix = '.so'
          if (!version.memory_bits eq 64) then suffix = '_64.so' $
          else                                 suffix = '.so'
        end
        else: suffix = '.so'
      endcase
    end
    'vms':     suffix = '.EXE'
    'Windows': suffix = '.DLL'
    'MacOS':   suffix = '.shlb'
    else: message, "Don't know what to do with: " + !version.os_family
  endcase
  libname = prefix + basename + suffix
end
