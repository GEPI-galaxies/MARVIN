	;-----------------------------------------------------------------------
	; FUNCTION GAMA_SMOOTH_CONTINUUM, Lambda,spectrum,mask_lines,pixels
	;
	; DESCRIPTION : Smoothing procedure using median filter 
	;				
	; 		INPUT : lambda, flux array and a mask vector ( 1-> masked pixel, 0-> unmasked )
	;				pixel -> size of the gaussian in pixels
	;		RETURN : Smooth spectrum vector
	; 	
	;----------------------------------------------------------------------
 
	FUNCTION GAMA_SMOOTH_CONTINUUM,Lambda,spectrum,mask,pixels
		lines=WHERE(mask EQ 1, COMPLEMENT=no_lines)
		spectrum_tmp=INTERPOL_V8(spectrum(no_lines),Lambda(no_lines),Lambda,/NAN)
		;savgolFilter = SAVGOL(150, 150, 0, 2)
		;Filter = DIGITAL_FILTER( 0, 0.085, 50, 5 )
		;continuum= CONVOL(spectrum_tmp, Filter, /EDGE_MIRROR)
		continuum=smooth(TS_smooth(spectrum_tmp,50 ),100 )
	
		return,continuum
	END
