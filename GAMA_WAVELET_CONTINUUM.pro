

	;-----------------------------------------------------------------------
	; FUNCTION GAMA_WAVELET_CONTINUUM, Lambda,spectrum,mask_lines,N_comp,order
	;
	; DESCRIPTION : Smoothing procedure using wavelet interpolation
	;				
	; 		INPUT : lambda, flux array and a mask vector ( 1-> masked pixel, 0-> unmasked )
	;				N_comp -> Number of wavelet orders to decompose
	;				order -> Number of low frequency order to merge	
	;		RETURN : Smooth spectrum vector
	; 	
	;----------------------------------------------------------------------
 
	FUNCTION GAMA_WAVELET_CONTINUUM,Lambda,spectrum,mask,N_comp,order,SILENT=SILENT
		lines=WHERE(mask EQ 1, COMPLEMENT=no_lines)
		
		spectrum_tmp=INTERPOL_V8(spectrum(no_lines),Lambda(no_lines),Lambda,/NAN)
		MR1D_ATROU,spectrum_tmp,decomposed,N_comp,/mirror
		order_max=N_comp-1
		order_min=N_comp-order
	
		tmp=decomposed[*,order_min:order_max]
		continuum=SUM(tmp,1)
		
			for j=0, 16 do begin
			spectrum_tmp(lines)=continuum(lines)
			MR1D_ATROU,spectrum_tmp,decomposed,N_comp;,/mirror
			tmp=decomposed[*,order_min:order_max]
			continuum=SUM(tmp,1)
			IF ~(KEYWORD_SET(SILENT)) then	print,".",format='(a,$)
		endfor
			IF ~(KEYWORD_SET(SILENT)) then begin 
			print," done",format='(a,$)
			print,""
			ENDIF
			
		return,continuum
	END


