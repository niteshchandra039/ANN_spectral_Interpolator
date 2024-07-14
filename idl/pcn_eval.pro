;+
; NAME:
;               PCN_EVAL
; PURPOSE: 
;               Evaluate a PCN interpolated model
;             
; USAGE:
;               interpolated = pcn_eval(pcn)
;
; ARGUMENTS:   
;   PCN:        PCN structure returned by the function PCN_INIT
;   PHYS:       Array of physical parameters
;
; KEYWORDS:
;   QUIET:      verbosity control
;
; DESCRIPTION: 
;   Return a structure containing the WCS and array of interpolated data
;
; HISTORY:
;               Creation Philippe Prugniel 2020/09/16
; 
;-
; CATEGORY:     PCN
;---------------------------------------------------------------------
function pcn_forward_propag, fun, input, weights

; I shall check if it is not faster to add a '1' in input and make a
; single matrix operation
  
  case fun of
     'linear' : begin
        output = weights[*,0] +  weights[*,1:*] # input
     end
     'logistic' : begin
        x = weights[*,0] +  weights[*,1:*] # input
        output = 1d / (1d + exp(-x))
     end
  endcase
    
  return, output
end

function pcn_eval, pcn, phys, QUIET=quiet

  if size(pcn, /TYPE) ne 8 then message, 'Argument PCN must be a structure'

  if n_elements(phys) ne n_elements(pcn.para) then message, 'Inconsistency beteen the number of parameters in the interpolator ('+n_elements(pcn.para)+') and the number of parameters in input ('+n_elements(phys)+')'
; shall check the presence of the tags expected in pcn
input = phys
  
; propagation through the layers
  case pcn.preproc of
     'identity': o = input  ; identity : input parameters = physical parameters
     'linear': o = input    ; I used to use linear instead of identity in the past
     'linscale': o = (input - pcn.preproc_para[*,0]) / pcn.preproc_para[*,1]
     else: message, 'preproc function not recognized ', + pcn.preproc
  endcase

  nlayer = n_elements(pcn.layer) ; number of layers

  for k=0,nlayer-1 do begin
     o = pcn_forward_propag(pcn.layer[k].afun, o, *pcn.layer[k].weights)
  endfor

  return, {crpix:pcn.crpix1, crval:pcn.crval1, cdelt:pcn.cdelt1, ctype:pcn.ctype1, data:o}
end
