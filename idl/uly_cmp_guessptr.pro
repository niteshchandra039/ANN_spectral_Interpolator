; the idea is to make a new general ulyss routine

function uly_cmp_guessptr, cmp

  if size(cmp, /TYPE) ne 8 then message, 'argument CMP is invalid, shall be a cmp structure'

  para = *cmp.para
  g = para[0].guess
  for k=1, n_elements(para)-1 do g = [g, para[k].guess]
  
; this may become part of a more general utilities to return
; information about a cmp, like the number of parameters, the list of their names, units, ...
  
; utility routines for cmp could be shared by all cmp and be useful to
; change limits, weight ... without directly manipulating the strucure
; (which is a bit tricky and maybe risky)
  
  return, g
end
