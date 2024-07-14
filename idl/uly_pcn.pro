;+
; NAME:
;                 ULY_PCN
; PURPOSE:
;                 Define a PCN component (perceptron interpolator)
; USAGE:
;                 cmp = uly_pcn(MODEL_FILE=model_file, LSF=lsf_file, NAME=name,      
;                               LIMITS=limits, FIXPAR=fixpar, WL=lim_weight)
;
; KEYWORDS:
;   MODEL_FILE    Name of the FITS file containing the PCN model.
;
;   LSF           Name of the file containing a relative LSF to be injected in 
;                 the template. 
;
;   NAME          A string, name of this component. By default a unique name is 
;                 generated.
;
;   LIMITS        [2,n] n is the number of free physical parameters. Each line of this 
;                 array is for a different physical parameter, and gives the low and 
;                 high bound of the range of valid values.
;                 By default these are the limits of the parameters, space defined in the file
;   
;   FIXPAR        [n] An array used to fix some the physical parameters,
;                 0 means that the parameter is free, and 1 that it is fixed.
;                 A fixed parameter keeps its initial value.
;                 By default all the parameters are free.
;          
;   WL            Limits for the weight of the component [min_w, max_w].
;                 The weight is in data units over component units.
;                 By default the weight is constrained to be positive.
;                 This constraint is applicable when multiple components are
;                 fitted (composite spectrum), it is ignored for a single
;                 component fit.
;
;   /REBIN_COEF   Computation mode, see the description below.
;       
; DESCRIPTION:
;     Return a structure describing a PCN component to be passed to program ULYSS. 
;
;     PCN stands for "perceptron" interpolator. It interpolates over a grid of models
;     (not necessarily a regular rectangular grid). It is described in Kumar et al. (20XX).
;     The interpolator is a multi-layer network fully described in the file.
;     
;     Some interpolators are available on the ULySS site.
;
;     ULY_PCN defines the component, i.e. it creates the structure cmp that will be later 
;     used to initialize and then evaluate a spectrum given a set of physical parameters.
;     See ULY_PCN_INIT about the initialization.
;
;
;     REBIN_COEF option:
;     Rather than rebinning the interpolated spectrum, it first rebin the 
;     coefficients, and then compute the interpolated spectrum. This option is not tested!
;     It is expected to be less precise, because of the interpolation of many 1D functions 
;     that are possibly very irregular. However it will be faster in the case when (i)
;     a series of several spectra are fitted with the same CMP, and (ii) when
;     the spectrum to analyse has a sampling much coarser than the interpolator.
;     (Note that the impact on the precision of the interpolated spectrum
;     is generally fully acceptable; it was the only available option before
;     the release 1.3 of ULySS.)
;
;
;     Controling the guesses
;     Properly choosing the values from which the minimization algorithm is started 
;     (the guesses) is important when the function has a complex behaviour. Randomly chosen 
;     guesses may result in stopping at a local minimum. 
;     By default the guesses are set at the middle of the parameter space.
;     This can be changed by calling the function ULY_CMP_GUESSPTR to return the array of 
;     pointers to the gueeses, and then alter the values. For example:
;            cmp = uly_pcn(file)
;            guess = uly_cmp_guessptr
;            *guess[1] = [1.0, 3.0]
;            ulyss, star, cmp, /PLOT
;     In this example, after defining the component, one retreive the guess pointer, and set 
;     guess of the second parameter to an array of two values. ULYSS then performs a global 
;     minimization (start from all the possible combinations of guesses, and retain the best solution).
;     
; OUTPUT:
;     cmp struct, for detailed explanation please check uly_fit.pro        
;
; REQUIRED FUNCTION:
;                 ULY_CMP_NAME
;                         
; EXAMPLE:
;     Fit a CFLIB star: 
;     cmp = uly_pcn(file)
;     star = uly_root+'/data/cflib_114642.fits'
;     ulyss, star, cmp, /PLOT
;
; HISTORY:
;                 Creation Ph. Prugniel 2022/01/09
;
;-
; CATEGORY:       ULY_PCN
;------------------------------------------------------------------------------
function uly_pcn, MODEL_FILE=model_file,                  $
                  LSF=lsf_file,                           $
                  NAME=name, LIMITS=limits,               $
                  FIXPAR=fixpar, WL=lim_weight,           $
                  REBIN_COEF=rebin_coef
                                                                    
compile_opt idl2
on_error, 2

common uly_path, uly_root

if n_elements(model_file) eq 0 then begin
    message, 'Keyword model_file must be a filename', /INFO
    return, 0
 endif

if size(model_file, /TYPE) ne 7 then begin
    message, 'Keyword model_file must be a filename', /INFO
    return, 0
endif
file = strtrim(model_file, 2)
if file_test(file) ne 1 then begin
    file += '.fits'
    if file_test(file) ne 1 then begin
        message, 'Keyword model_file must be a filename ('+model_file+')', /INFO
        return, 0
    endif
endif

if n_elements(rebin_coef) eq 0 then rebin_coef = 0

s_wl = size(lim_weight)
if s_wl[0] gt 0 then $
   if not s_wl[0] eq 2 then $
   message, 'The weight limits have to be of the type arr(2)'

; create the struct describing PCN component
init_data = pcn_init(file)

descr = ''
if n_elements(file) eq 1 then descr += 'model:' + file + ' ' 
if n_elements(lsf_file) eq 1 then descr += 'lsf:' + lsf_file + ' '
if rebin_coef eq 1 then descr += 'rebin_coef ' 

namec = uly_cmp_name(name)

cmp = {name:namec, $
       descr:descr, $
       init_fun:'uly_pcn_init', $
       init_data:ptr_new(init_data), $
       eval_fun:'', $
       eval_data:ptr_new(), $
       para:ptr_new(/ALLO), $
       start:0d, $
       step: 0d, $
       npix: 0l, $
       sampling:-1s, $
       mask:ptr_new(), $
       weight:0d, $
       e_weight:0d, $
       l_weight:0d, $
       lim_weig:(machar(/DOUBLE)).xmax*[0,1] $
      }


if n_elements(lim_weight) gt 0 then cmp.lim_weig = double(lim_weight)

; load the parameters and data in the component struct
npara = size(init_data.para, /N_ELEMENTS) ; number of parameters
*cmp.para = replicate({name:'', unit:'', guess:ptr_new(), step:1D-2, $
                       limits:[0d,0d], limited:[1,1], fixed:0S, $
                       value:0D, error:0D, dispf:''}, $
                      npara)       

(*cmp.para).name = (init_data.para).name
(*cmp.para).unit = (init_data.para).unit
(*cmp.para).limits = (init_data.para).limits

if n_elements(limits) ne 0 then begin
   if size(limits, /N_DIM) ne 2 then message, 'LIMITS keyword has invalid dimensions, it is expected to be a 2D array'
   sz = size(limits, /DIM)
   if sz[0] ne 2 or sz[1] ne npara then message, 'LIMITS keyword has invalid dimension'
   (*cmp.para).limits = limits
endif

; Set default guess (average of lower and upper limits)
for k=0,npara-1 do (*cmp.para)[k].guess = ptr_new(((*cmp.para)[k].limits[0]+(*cmp.para)[k].limits[1])/2d)

; Set the fixed parameters
if n_elements(fixpar) gt 0 then (*cmp.para)[0:n_elements(fixpar)-1].fixed = fixpar

; Set the steps (for derivation)
(*cmp.para).step = (init_data.para).step  

; Notes:
; - Still missing "dispf" in pcn
;(*cmp.para)[0].dispf = 'exp'  ; says to display exp(para)

return, cmp

end
