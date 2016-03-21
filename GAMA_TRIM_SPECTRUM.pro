
	;-----------------------------------------------------------------------
	; FUNCTION TRIM_SPECTRUM
	; DESCRIPTION : Trim the spectra to a given boundary 
	;               Use sigma cleaping to calculate the mean flux.
	;---------------------------------------------------------------------
	
	FUNCTION GAMA_TRIM_SPECTRUM,lambda,spectra,region2trim
	spectrum_trim=spectra(WHERE(lambda GE region2trim[0] and lambda LE region2trim[1]))
	return, spectrum_trim
	END
