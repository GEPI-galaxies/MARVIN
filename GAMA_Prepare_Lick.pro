FUNCTION GAMA_Prepare_Lick, struct_spec, TRIM=TRIM, EXTINCTION=EXTINCTION, MASK_SLICER=MASK_SLICER, MASK_EMI=MASK_EMI

@GAMA_Main_Config

lambda=struct_spec.lambda
spectrum_orig=struct_spec.Spectrum
error_orig=struct_spec.Error
sky_orig=struct_spec.Sky
objet=struct_spec.name


IF TAG_EXIST(struct_spec,'Sensibility') then begin
Sensibility_orig=struct_spec.Sensibility
ENDIF ELSE begin
Sensibility_orig=fltarr(n_elements(lambda))
Sensibility_orig(*)=1.
ENDELSE

sxaddpar,struct_spec.header,'HISTORY','Process by MARVIN'

IF KEYWORD_SET(TRIM) then begin
	;TRIM if too noisy spectrum at the beginning of the spectrum 
		max_value_lamb=max(lambda)- 100.
		min_value_lamb=min(lambda)+ 100.
		spectrum_orig=GAMA_TRIM_SPECTRUM(lambda,spectrum_orig,[min_value_lamb,max_value_lamb])
		error_orig=GAMA_TRIM_SPECTRUM(lambda,error_orig,[min_value_lamb,max_value_lamb])
		sky_orig=GAMA_TRIM_SPECTRUM(lambda,sky_orig,[min_value_lamb,max_value_lamb])
		Sensibility_orig=GAMA_TRIM_SPECTRUM(lambda,Sensibility_orig,[min_value_lamb,max_value_lamb])
		lambda=GAMA_TRIM_SPECTRUM(lambda,lambda,[min_value_lamb,max_value_lamb])
		sxaddpar,struct_spec.header,'HISTORY','Trim bellow '+string(min_value_lamb)
		;ENDIF	
ENDIF


IF KEYWORD_SET(EXTINCTION) then begin
	;Correct from galatic extinction
	;Read value of the extinction from schlegel extinction map and correct from foreground dust
	RA=SXPAR(struct_spec.header,'RA')
	DEC=SXPAR(struct_spec.header,'DEC')
	euler, RA, DEC, l, b, 1 ;Convert to galatic coordinates
	EBV = dust_getval(l, b,ipath=PATH_DUST,/interp)
	CCM_UNRED, lambda, spectrum_orig, EBV
	sxaddpar,struct_spec.header,'EBV',EBV	
ENDIF	
	redshift=SXPAR(struct_spec.header,'Z')
	NAXIS1=SXPAR(struct_spec.header,'NAXIS1')
	lambda_o_rest=lambda/(1+redshift)
	spectrum_orig=spectrum_orig*(1+redshift)

	Dim_spectrum = n_elements(spectrum_orig)
	;Create mask ( 10 -> NAN pixels, 6 -> Slicer region)
	mask=intarr(Dim_spectrum)
	mask(*)=0
	nan_index = WHERENAN(spectrum_orig,n_nan) ;Mask NAN regions
		IF n_nan GT 1 then begin 
		spectrum_orig(nan_index)=0.
		error(nan_index)=0.
		mask(nan_index)=10
		ENDIF

	IF KEYWORD_SET(MASK_SLICER) then begin
	slicer=WHERE(lambda_o_rest GT 5500.*(1+redshift) and lambda_o_rest LT 6000.*(1+redshift),cnt_slicer ) ;Mask slicer region
	IF cnt_slicer GT 1 then mask(slicer)=6
	ENDIF


new_struct_spec={Lambda:lambda_o_rest,Spectrum:spectrum_orig,Error:error_orig,Sensibility:Sensibility_orig,Sky:sky_orig,Mask:mask,name:struct_spec.name, header:struct_spec.header}


return, new_struct_spec


END