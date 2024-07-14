;+
; NAME:
;               PCN_INIT
; PURPOSE: 
;               Read a PCN file
;             
; USAGE:
;               pcn = pcn_init(filename [, /QUIET])
;
; ARGUMENTS:   
;   FILENAME:   Name of the PCN file (FITS format)
;
; KEYWORDS:
;   QUIET:      verbosity control
;
; DESCRIPTION: 
;      Read a PCN FITS file and load the coefficients and meta-data into
;      a pcn structure.
;      The PCN format is specified in Kumar et al. 20XX.
;
; HISTORY:
;               Creation Philippe Prugniel 2020/09/16
; 
;-
; CATEGORY:     PCN
;---------------------------------------------------------------------
function pcn_init, filename, QUIET=quiet

compile_opt idl2, hidden
on_error, 2

if n_elements(filename) eq 0 then begin
    message, 'no input model file name...', /INFO
    print, 'Usage: pcn_init, <filename> [, /QUIET]'
    return, 0
endif

fits_read, filename, wg, hdr, EXTEN_NO=0, MESSAGE=mess0  

; The primary HDU contains meta-data describing the interpolator, and
; the terminal weights, linking the last hidden layer to the output layer

; Check if the file is PCN in a recognised version
pcn_version = strtrim(sxpar(hdr, 'I_VERSIO', /SILENT, COUNT=cnt), 2)
if cnt eq 0 then begin
   message, 'This is apparently not a PCN file (miss kw I_VERSIO)'
endif
if strmid(pcn_version, 0, 2) ne '0.' then begin
   message, 'The PCN version is not recognized (I_VERSIO='+pcn_version+')'
endif

; Number of hidden layers
pcn_hlayer = sxpar(hdr, 'I_HLAYER', /SILENT, COUNT=cnt)

if cnt eq 0 then begin
   message, 'Error in PCN file (kw I_HLAYER missing)'
endif

; pre-processing function
; I_PREPRO is the name of the pre-processing function (transforming
; the physical parameters into input parameters), by default we
; assume an identity (input parameters = physical parameters)
pcn_preproc = strlowcase(strtrim(sxpar(hdr, 'I_PREPRO', /SILENT, COUNT=cnt),2))
if cnt eq 0 then pcn_preproc = 'identity' $
else if strlen(pcn_preproc) eq 0 then pcn_preproc = 'identity' 

case pcn_preproc of
   'identity': preproc_para = 1
   'linscale': begin
      fits_read, filename, data, h, EXTEN_NO=pcn_hlayer, MESSAGE=mess, /HEADER_ONLY
      npara = sxpar(h, 'NAXIS2', /SILENT, COUNT=cnt)
      if cnt ne 1 then message, 'no NAXIS2 in last extension?'
      npara -= 1 ; number of physical parameters
      ; this is actually the number of parameters, which in the case of linscale appears to be the number of PHYSICAL parameters
      
; Note that we shall have a desciption of the physical parameters in
; the PHDU header, so, we shall get npara from that (no need to read the last extension)
      
      ; currently hardcoding to the case of the toy interpolator

; Read the bias (I_PR_B) and scale (I_PR_S) for each param
; something more flexible (for the next version)
      preproc_coef_b = sxpar(hdr, 'I_PR_B*', COUNT=cnt)
      if cnt ne npara then message, 'Incorrect number of I_PR_Bnn keywords, expect ' + npara + 'got ' + cnt
      preproc_coef_s = sxpar(hdr, 'I_PR_S*', COUNT=cnt)
      if cnt ne npara then message, 'Incorrect number of I_PR_Snn keywords, expect ' + npara + 'got ' + cnt
      preproc_para = [[preproc_coef_b], [preproc_coef_s]]            
;      preproc_para = dblarr(npara,2)
;      preproc_para[0,0] = double(string(sxpar(hdr, 'I_PR_TMU')))
;      preproc_para[0,1] = double(string(sxpar(hdr, 'I_PR_TSI')))
;      preproc_para[1,0] = double(string(sxpar(hdr, 'I_PR_GMU')))
;      preproc_para[1,1] = double(string(sxpar(hdr, 'I_PR_GSI')))
;      preproc_para[2,0] = double(string(sxpar(hdr, 'I_PR_MMU')))
;      preproc_para[2,1] = double(string(sxpar(hdr, 'I_PR_MSI')))      
      end
endcase
   
; post-processing function
pcn_postproc = strlowcase(strtrim(sxpar(hdr, 'I_POSTPR', /SILENT, COUNT=cnt),2))
if cnt eq 0 then pcn_postproc = 'identity' $
else if strlen(pcn_postproc) eq 0 then pcn_postproc = 'identity' 

; Read the WCS
naxis1 = sxpar(hdr, 'NAXIS1', /SILENT, COUNT=cnt)
crpix1 = sxpar(hdr, 'CRPIX1', /SILENT, COUNT=cnt)
if cnt eq 0 then begin
   message, /INFO, 'Missing CRPIX1 keyword, assume CRPIX1=1'
   crpix1 = 1d
endif else crpix1 = double(string(crpix1)) 
crval1 = double(string(sxpar(hdr, 'CRVAL1')))
cdelt1 = double(string(sxpar(hdr, 'CDELT1')))

ctype1   = strtrim(sxpar(hdr, 'CTYPE1', /SILENT, COUNT=cnt), 2)
if cnt eq 0 then ctype1 = 'AWAV'

; Load the physical parameters meta data

; Physical parameters
; I_P_NM : name
; I_P_UN : unit
; I_P_LB : low bound
; I_P_HB : high bound
; I_P_ST : step for numerical derivation
name = string(sxpar(hdr, 'I_P_NM*', COUNT=npara))
unit = string(sxpar(hdr, 'I_P_UN*', COUNT=cnt))
lb = string(sxpar(hdr, 'I_P_LB*', COUNT=cnt))
hb = string(sxpar(hdr, 'I_P_HB*', COUNT=cnt))

stp = double(string(sxpar(hdr, 'I_P_ST*', COUNT=cnt))) ; derivation steps

; TEMPORARY: The header does not yet contain I_P_NM*, .... so we are HARDCODING
; The following information will have to be read from keywords
if npara eq 0 then begin
   message, /INFO, 'Old data model, kws I_P_NM, I_P_UN, I_P_LB, I_P_HB, I_P_ST missing)'
   npara = 3
   name = ['Teff', 'Logg', '[Fe/H]']
   unit = ['[K]', 'cm/s2', 'dex']
   lb = [3500., 0., -3.]
   hb = [12000., 6., 0]
   stp = [1d, 1d-2, 1d-2]
endif

;para[0].dispf = 'exp'  ; says to display exp(para)

para = replicate ({name:'', unit:'', $
                   limits:[0d,0d], step:0d}, $
                  npara) 

para.name = name
para.unit = unit
para.limits = [transpose(lb),transpose(hb)] 
para.step = stp

; Create the layer array to contain the weights and metadata
layer = replicate({weights:ptr_new(/ALLO), afun:''}, pcn_hlayer+1)

; Load the weights and metadata
for k=0,pcn_hlayer do begin   

   ; Index of the target layer
   pcn_ilayer = sxpar(hdr, 'I_LAYER', /SILENT, COUNT=cnt)
   if cnt eq 0 then begin
      message, 'Error in PCN file (kw I_LAYER missing)'
   endif

   ; Name of the activation function
   pcn_afun = strtrim(sxpar(hdr, 'I_AFUNC', /SILENT, COUNT=cnt),2)
   if cnt eq 0 then begin
      message, 'Error in PCN file (kw I_AFUNC missing)'
   endif
   layer[pcn_ilayer-1].afun = pcn_afun

   ; Attach the weights
   layer[pcn_ilayer-1].weights = ptr_new(wg)
   if k lt pcn_hlayer then fits_read, filename, wg, hdr, EXTEN_NO=k+1, MESSAGE=mess
endfor

return, {preproc:pcn_preproc, preproc_para:preproc_para, postproc:pcn_postproc, layer:layer, naxis1:naxis1, crpix1:crpix1, crval1:crval1, cdelt1:cdelt1, ctype1:ctype1, para:para}

end

