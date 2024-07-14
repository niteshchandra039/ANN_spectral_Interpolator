pro param


spec_files = FILE_SEARCH('spectra/*.txt')


model_path = '../trial_6/toy_interpolator_v0.0.3_trial_6.fits'


cmp = uly_pcn(MODEL = model_path)
guess = uly_cmp_guessptr(cmp)

for i=0, 5 do begin
    print, i
	spec_file = spec_files[i]
	wl_file = 'wl/'+ strmid(spec_file, 8, STRLEN(spec_file))
	id = strmid(spec_file, 49, 9)
	spec_data = read_ascii(spec_file)
	spec_data = spec_data.(0)

	wl_data = read_ascii(wl_file)
	wl_data = wl_data.(0)

	spec = uly_spect_alloc(TITLE=spec_file, data=spec_data, sampling=0, start=wl_data[0], step=wl_data[1]-wl_data[0]) 


	; *guess[0] = strmid(star_name, 4,4)
	; *guess[1] = strmid(star_name, 9,4)
	; *guess[2] = strmid(star_name, 13,4)

	ulyss, spec, cmp, FILE_OUT='output_temp/'+id, /PLOT, /QUIET
	uly_spect_free, spec
	heap_gc, /PTR
endfor

end
