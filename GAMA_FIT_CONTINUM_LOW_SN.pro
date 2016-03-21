
;-----------------------------------------------------------------------
; PROCEDURE GAMA_FIT_CONTINUM_LOW_SN
; DESCRIPTION : Fit the stellar continuum with wavelet 
;               for low S_N spectrum
;---------------------------------------------------------------------

FUNCTION GAMA_FIT_CONTINUM_LOW_SN,lambda_rest,spectrum

@GAMA_Main_Config

	
	;Create mask ( 10 -> NAN pixels, 6 -> Slicer region)
	Dim_spectrum=n_elements(lambda_rest)
	mask=intarr(Dim_spectrum)
	mask(*)=0
	nan_index = WHERENAN(spectrum,n_nan) ;Mask NAN regions
		IF n_nan GT 1 then begin 
		spectrum(nan_index)=0.
		;error(nan_index)=0.
		mask(nan_index)=1
		ENDIF

		; Mask spikes
		Mask_Spikes=GAMA_WAVELET_SPIKES(spectrum,9,2,3)
		
		; Mask emission lines 
		Mask_lines=GAMA_LOAD_emiMask(lambda_rest)
		all_mask=Mask_lines+mask+Mask_Spikes
		tmp=WHERE(all_mask GT 1,cnt)
			IF cnt GT 1 then all_mask(tmp)=1
		
		continuum=GAMA_WAVELET_CONTINUUM(lambda_rest,spectrum,all_mask,11,3,SILENT=SILENT)
		;oplotObj1 = Obj_New('cgOverPlot',lambda_rest,continuum,thick=2, color="red") 
		;oplotObj2 = Obj_New('cgOverPlot',lambda_rest,all_mask,thick=1, color="green") 
	;	oplotObj3 = Obj_New('cgOverPlot',lambda_rest,Spikes,thick=1, color="pink") 
		
		;cgplot,lambda_rest,smooth(spectrum,3), OPLOTS=[oplotObj1,oplotObj2],xrange=[3700.,6700.]
	
		 return,continuum
END

