;+
; NAME:            ULY_PCN_EVAL
; 
; PURPOSE:           
;                  Evaluate a PCN model spectrum.
; 
; USAGE:
;                  array = uly_tgm_eval(eval_data, para)
; 
; ARGUMENTS:
;   EVAL_DATA:     Input, structure describing the interpolator
;
;   PARA:          Input, array of physical parameters
;
; DESCRIPTION:
;      Evaluate a PCN interpolated spectrum.
;
;      This function is automatically called by the minimization routines:
;      ULY_FIT, ULY_FIT_LIN to evaluate a PCN model (see ULY_PCN).
; 
; OUTPUT:          
;                  Evaluated stellar model spectrum, given as an array
;
; AUTHOR:
;                  Philippe Prugniel, 2020/09/24
;
;-
; CATEGORY:        ULY_PCN
;------------------------------------------------------------------------------
function uly_pcn_eval, eval_data, para

sp = pcn_eval(*eval_data.pcn, para)

pcn_spectrum = sp.data

; Resample in wavelength if necessary (when option REBIN_COEF is not given)
if eval_data.sampling ne eval_data.mod_samp or $
   eval_data.start ne eval_data.mod_start or $
   eval_data.step ne eval_data.mod_step or $
   eval_data.npix ne n_elements(pcn_spectrum) then begin
   spec = uly_spect_alloc(DATA=pcn_spectrum, START=eval_data.mod_start, $
                          STEP=eval_data.mod_step, SAMPLING=eval_data.mod_samp)
   wrange = eval_data.start + [0d, eval_data.npix * eval_data.step]
   if eval_data.sampling eq 1 then wrange = exp(wrange)
   c = 299792.458d              ; Speed of light in km/s
   velscale = eval_data.step * c
   if eval_data.sampling eq 1 then $
      spec = uly_spect_logrebin(spec, velscale, WAVERANGE=wrange, /OVER) $
   else $
      spec = uly_spect_linrebin(spec, eval_data.step, WAVERANGE=wrange, /OVER)
   pcn_spectrum = *spec.data
   uly_spect_free, spec
endif

;; Convolve LOSVD in case giving lsf_file
;; (if REBIN_COEF is given, the LSF injection shall be made at initialization.)
;if n_elements(eval_data.lsf) gt 0 and eval_data.lsf ne 'no_lsf' then begin
;    spec = uly_spect_alloc(DATA=tgm_model_evalhc, START=eval_data.start, $
;                           STEP=eval_data.step, SAMPLING=1)
;    uly_spect_lsfconvol, eval_data.lsf, spec, /QUIET
;    tgm_model_evalhc = *spec.data
;    uly_spect_free, spec
;endif

su = uly_spect_alloc(DATA=pcn_spectrum, START=eval_data.start, STEP=eval_data.step, SAMPLING=eval_data.sampling)

;uly_spect_plot, su

return, pcn_spectrum

end

