pro param_new


filename = FILE_SEARCH('*.txt')

;for j=0, 7 do begin
j=8
	trial='trial_'+strtrim(j, 1)
	path = '/home/nitesh/nitesh/PhD/GSL_Interpolator/Trials_for_different_architecture/'+trial+'/'
	file = path+'output_trial/temp_'+trial+'_3.dat'
	cmp = uly_pcn(MODEL=path+'toy_interpolator_v0.0.3_'+trial+'.fits')
	guess = uly_cmp_guessptr(cmp)

	for i= 0, 2907 do begin
		print, trial, i
		star = filename(i)
		data = read_ascii(star)
		data = data.(0)
		wl = read_ascii('/home/nitesh/nitesh/spectral_library/GSL/pheonix_txt/wl/'+strmid(star,0, 64)+'wl')
		wl = wl.(0)
		spec = uly_spect_alloc(TITLE=star, data=data, sampling=0, start=wl[0], step=wl[1]-wl[0])
		*guess[0] = strmid(star, 4,4)
		*guess[1] = strmid(star, 9,4)
		*guess[2] = strmid(star, 13,4)
		ulyss, spec, cmp, FILE_OUT=path+'output_trial/pheonix_'+trial+'_3', MD=0, KM=0, /QUIET
		uly_spect_free, spec
		heap_gc, /PTR
		openw, 21, file, /append
		printf, 21, i,',',star
		close, 21
	;endfor
endfor
heap_gc, /PTR
end

