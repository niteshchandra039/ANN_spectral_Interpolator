;+
; NAME:
;               ULY_PCN_INIT
; PURPOSE: 
;               Initialize a PCN component (perceptron)
;             
; USAGE:
;               cmp_new = uly_pcn_init (cmp, WAVERANGE=lamrange, VELSCALE=velscale, QUIET=quiet)
;
; ARGUMENTS:   
;   CMP:        Component defined by using ULY_PCN.
;
; KEYWORDS:
;   WAVERANGE:  Wavelength range used when the model will be
;               evaluated with ULY_PCN_EVAL.
;
;   VELSCALE:   Size of one pixel in km/s. Used when the model will be
;               evaluated and log-rebinned with ULY_PCN_EVAL.
;
;   QUIET:      verbosity control
;
; DESCRIPTION: 
;      This initialization function is automatically executed by 
;      ULY_FIT_INIT when a component has been defined by using ULY_PCN.
;
;      The name of this function is recorded as a member of
;      cmp.init_fun in the TGM component structure (see ULY_FIT for a
;      complete description of the component structure).
;
;      ULY_PCN_INIT loads the coefficients from a FITS file whose name is saved in
;      *cmp.init_data, the WCS
;
;
; HISTORY:
;               Creation Philippe PRUGNIEL 2020/09/24
; 
;-
; CATEGORY:     ULY_PCN
;---------------------------------------------------------------------

function uly_pcn_init, cmp, WAVERANGE=lamrange, VELSCALE=velscale, $
                            SAMPLING=sampling, STEP=step, QUIET=quiet
compile_opt idl2
on_error, 2

init_data = *cmp.init_data

cmp.eval_fun = 'uly_pcn_eval'

;; In the present version the coeff arrays are yet in 'init_data' (loaded whe uly_pcn is called)
;; This may not be optimum as it implies that the interpolator's data are stored twice, once in init_data, and second in eval_data

;if size(init_data.model, /TYPE) eq 7 then begin
;   if file_test(init_data.model) ne 1 then begin
;      print, 'Error, model file does not exist!'
;   endif
;endif

naxis1 = init_data.naxis1
ctype1 = init_data.ctype1
crpix1 = init_data.crpix1
crval1 = init_data.crval1
cdelt1 = init_data.cdelt1

if ctype1 eq 'AWAV' then pcn_samp = 0 $
else if ctype1 eq 'AWAV-LOG' then pcn_samp = 1 $
else message, 'Unsupported CTYPE1 value ' + ctype1

if ptr_valid(cmp.eval_data) then begin
;  this cmp was yet initialized, we shall check whether we have to renitialize
;    compare actual arguments with those stored in eval_data
   arg = (*cmp.eval_data).arg
   init = 1
   if n_elements(lamrange) gt 0 then begin
      if n_elements(lamrange) ne n_elements(arg.lamrange) then init *= 0 else $
         if total(lamrange eq arg.lamrange) ne n_elements(lamrange) then init *= 0
   endif else if arg.lamrange ne -1 then init *= 0
   if n_elements(velscale) eq 1 then if velscale ne arg.velscale then init *= 0 else if arg.velscale ne -1 then init *= 0
   if n_elements(sampling) eq 1 then if sampling ne arg.sampling then init *= 0 else if arg.sampling ne -1 then init *= 0
   if n_elements(step) eq 1 then if step ne arg.step then init *= 0 else if arg.step ne -1 then init *= 0

   if init eq 1 then return, cmp ; no need to re-initialize
endif


; attach the terminal weights in a spect structure that we will trim and/or resample
nlayer = n_elements(init_data.layer)
s = uly_spect_alloc(START=crval1, STEP=cdelt1, SAMPLING=pcn_samp)
s.data = init_data.layer[nlayer-1].weights  ; terminal weights

;; attention bug: we attach the pointer from init_data, and we modify the data
;; therefore if the model's wavelength range was trimmed, we cannot 're-initialize'

wr = crval1 + [0, (naxis1-1)*cdelt1]
if ctype1 eq 'AWAV-LOG' then wr = exp(wr)
if n_elements(lamrange) gt 0 then wr[0] = max([lamrange[0],wr[0]])
if n_elements(lamrange) gt 1 then wr[1] = min([lamrange[1],wr[1]])

; rebin_coef is currently not loaded in init_data (init_data is directly given by the pcn routine, and rebin_coef is a ULySS thing)
;if init_data.rebin_coef eq 1 then begin ; resample the terminal weights
; for the moment we NEVER REBIN
if 0 eq 1 then begin            ; resample the terminal weights
   sampl = 1
   if n_elements(sampling) eq 1 then sampl = sampling 
   if sampl eq 1 then begin
      undefine, velsc
      if n_elements(velscale) eq 1 then velsc = velscale $
      else if n_elements(step) eq 1 then velsc = step * 299792.458d 
      s = uly_spect_logrebin(s, velsc, WAVERANGE=wr, /OVER) 
   endif else if sampl eq 0 then begin
      s = uly_spect_linrebin(s, step, WAVERANGE=wr, /OVER) 
   endif
endif else begin
   s = uly_spect_extract(s, WAVERANGE=wr, /OVER)
endelse


; Description of the model coefs stored in cmp (wavelength range, sampling and step)
mod_samp = s.sampling
mod_start = s.start
mod_step = s.step

; Determine the WCS of the cmp itself
cmp.sampling = 1 ; by default the cmp will be log-sampled
if n_elements(sampling) eq 1 then cmp.sampling = sampling

if cmp.sampling eq mod_samp then cmp.step = mod_step 
if n_elements(velscale) eq 1 then begin
   if cmp.sampling ne 1 then message, 'Inconsistency in the arguments'
   cmp.step = velscale/299792.458d0
endif
if n_elements(step) eq 1 then cmp.step = step

resam = 0
if n_elements(lamrange) gt 0 then if wr[0] ne lamrange[0] then resam = 1
if n_elements(lamrange) gt 1 then if wr[1] ne lamrange[1] then resam = 1

if mod_samp ne cmp.sampling or mod_step ne cmp.step or resam eq 1 then begin
   s1d = uly_spect_extract(s, ONED=0) ; to compute the WCS
   if cmp.sampling eq 0 then begin
      undefine, step
      if cmp.step ne 0 then step= cmp.step
      s1d = uly_spect_linrebin(s1d, step, WAVERANGE=lamrange, /OVER)
   endif else if cmp.sampling eq 1 then begin
      undefine, velscale
      if cmp.step ne 0 then velscale = cmp.step * 299792.458d0
      s1d = uly_spect_logrebin(s1d, velscale, WAVERANGE=lamrange, /OVER)
   endif else message, 'Cannot yet resample to sampling=2'
   cmp.start = s1d.start
   cmp.step = s1d.step
   cmp.npix = (size(*s1d.data, /DIM))[0]
   cmp.sampling = s1d.sampling
   heap_free, s1d
endif else begin
   cmp.start = s.start
   cmp.step = s.step
   cmp.npix = (size(*s.data, /DIM))[0]
   cmp.sampling = s.sampling
endelse

; We currently do not have a "goodpixel" list in PCN
;if n_elements(*s.goodpix) gt 0 then begin
;   ptr_free, cmp.mask
;   cmp.mask = ptr_new(bytarr(cmp.npix)) 
;   (*cmp.mask)[*s.goodpix] = 1 
;endif

s.data = ptr_new() 
heap_free, s

ptr_free, cmp.eval_data

;if n_elements(*init_data.lsf_file) gt 0 then lsf = *init_data.lsf_file $
;else
lsf = 'no_lsf'

; store the arguments in eval_data (memorize how init was made)
if n_elements(lamrange) eq 1 then arg_lamrange = lamrange else arg_lamrange = -1
if n_elements(velscale) eq 1 then arg_velscale = velscale else arg_velscale = -1
if n_elements(sampling) eq 1 then arg_sampling = sampling else arg_sampling = -1
if n_elements(step) eq 1 then arg_step = step else arg_step = -1
arg = {lamrange:arg_lamrange, velscale:arg_velscale, sampling:arg_sampling, step:arg_step}

cmp.eval_data = ptr_new({pcn:cmp.init_data, $
                        start:cmp.start, step:cmp.step, npix:cmp.npix, sampling:cmp.sampling, $
                        mod_samp:mod_samp, mod_start:mod_start, mod_step:mod_step, $
                        lsf:lsf, arg:arg})

return, cmp

end

