
	;-----------------------------------------------------------------------
	; FUNCTION GAMA_SLIDING_MOMENTS, spectrum , win_width
	;
	; DESCRIPTION : Calculate the sliding mean and sigma of a signal 
	;				using a box window of n_window size
	; 		INPUT : spectrum
	;				win_width - Sliding window half width
	; 				slide_incr - Slide for each iteration		
	;		RETURN : mean and sigma spectrum 
	; 	
	;----------------------------------------------------------------------
 
	FUNCTION GAMA_SLIDING_MOMENTS,spectrum,win_width
		
		Mean_array=fltarr(n_elements(spectrum))
		Sigma_array=fltarr(n_elements(spectrum))
		slide_incr=1
		for i=0, n_elements(spectrum) -1 do begin
		
			IF i LT win_width then begin
			;left_mirror=REVERSE(spectrum( 0 : win_width -i))
			;n_mirror=n_elements(left_mirror)
			;spectrum=[left_mirror,spectrum]
			;i_new=n_mirror+i
			;Mean_array[i]=mean(spectrum(i_new-win_width:i_new+win_width))
			;Sigma_array[i]=sigma(spectrum(i_new-win_width:i_new+win_width))	
			Mean_array[i]=mean(spectrum(0:i+win_width))
			Sigma_array[i]=sigma(spectrum(0:i+win_width))	

			endif 
		    IF i GE n_elements(spectrum) - win_width -1 then begin
		   	right_mirror=REVERSE(spectrum( n_elements(spectrum) -1 - win_width : n_elements(spectrum) -1))
			n_mirror=n_elements(right_mirror)
			spectrum=[spectrum,right_mirror]
			Mean_array[i]=mean(spectrum(i-win_width:i+win_width))
			Sigma_array[i]=sigma(spectrum(i-win_width:i+win_width))	
		    endif 
		    IF i GE win_width and i LT n_elements(spectrum) - win_width -1 then begin
			Mean_array[i]=median(spectrum(i-win_width:i+win_width))
			Sigma_array[i]=sigma(spectrum(i-win_width:i+win_width))
			endif
		endfor	
		return,[[Mean_array],[Sigma_array]]
	END
