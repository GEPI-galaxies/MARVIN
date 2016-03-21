
	
	;-----------------------------------------------------------------------
	; FUNCTION GAMA_WAVELET_SPIKES,spectrum,N_comp,order
	;
	; DESCRIPTION : Smoothing procedure using wavelet interpolation
	;				
	; 		INPUT : lambda, flux array and a mask vector ( 1-> masked pixel, 0-> unmasked )
	;				N_comp -> Number of wavelet orders to decompose
	;				order -> Number of low frequency order to merge	
	;		RETURN : Smooth spectrum vector
	; 	
	;----------------------------------------------------------------------
 
	FUNCTION GAMA_WAVELET_SPIKES,spectrum,N_comp,order,CLIPSIG
		
		Mask_spikes=fltarr(n_elements(spectrum))
		
		Mask_spikes(*)=0.
		MR1D_ATROU,spectrum,decomposed,N_comp,/mirror
		order_max=order
		order_min=0
	
		tmp=decomposed[*,order_min:order_max]
		spikes=SUM(tmp,1)
		spikes=spectrum
		Despike=spectrum-Spikes
		
		moments_array=GAMA_SLIDING_MOMENTS(spikes,20)
		mean_array=smooth(moments_array[*,0],20)
		sigma_array=smooth(moments_array[*,1],20)
		spikes_position = WHERE( abs(Spikes - mean_array) GT  sigma_array *CLIPSIG)
		;MEANCLIP, Despike, Mean_rob, Sig_rob,CLIPSIG=CLIPSIG,SUBS=no_clip
		;spikes_position=mg_complement(no_clip, n_elements(spectrum))
		;IF n_elements(spikes_position) GT 1 then Mask_spikes(spikes_position)=1.
		;cgplot, Spikes - mean_array
		;cgoplot,mean_array, color='red'
		;cgoplot,CLIPSIG*sigma_array, color='grey'
		;cgoplot,Mask_spikes*max(sigma_array),color='pink',thick=3

		return,Mask_spikes
	END
