FUNCTION GAMA_PREPARE_STARLIGHT, struct_spec, TRIM=TRIM, EXTINCTION=EXTINCTION, MASK_SLICER=MASK_SLICER, MASK_EMI=MASK_EMI , MASK_BLUE=MASK_BLUE

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

print, "   [STEP 0] : Preparing spectrum"
sxaddpar,struct_spec.header,'HISTORY','Process by MARVIN'

IF KEYWORD_SET(TRIM) then begin
	;TRIM if too noisy spectrum at the beginning of the spectrum 
	;relative_error=abs(error_orig/spectrum_orig)
	;mean_error_bleu=mean(relative_error(WHERE( lambda LT 4000)),/NAN)
	;	IF mean_error_bleu GT 3. then begin
		max_value_lamb=max(lambda)- 100.
		min_value_lamb=min(lambda)+ 100.
		spectrum_orig=GAMA_TRIM_SPECTRUM(lambda,spectrum_orig,[min_value_lamb,max_value_lamb])
		error_orig=GAMA_TRIM_SPECTRUM(lambda,error_orig,[min_value_lamb,max_value_lamb])
		sky_orig=GAMA_TRIM_SPECTRUM(lambda,sky_orig,[min_value_lamb,max_value_lamb])
		Sensibility_orig=GAMA_TRIM_SPECTRUM(lambda,Sensibility_orig,[min_value_lamb,max_value_lamb])
		lambda=GAMA_TRIM_SPECTRUM(lambda,lambda,[min_value_lamb,max_value_lamb])
		sxaddpar,struct_spec.header,'HISTORY','Trim bellow '+string(min_value_lamb)
	print, "        -Trim spectrum under 4000ang: done"
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

	;Normalize to Starlight sampling (1A)
	lambda_min=floor(min(lambda_o_rest, /NAN))
	lambda_max=ceil(max(lambda_o_rest,/NAN))

	lambda_rest=indgen(lambda_max-lambda_min)+lambda_min
	spectrum=INTERPOL_V8(spectrum_orig,lambda_o_rest,lambda_rest,/NAN)
	error=INTERPOL_V8(error_orig,lambda_o_rest,lambda_rest,/NAN)
	sky=INTERPOL_V8(sky_orig,lambda_o_rest,lambda_rest,/NAN)
	Sensibility=INTERPOL_V8(Sensibility_orig,lambda_o_rest,lambda_rest,/NAN)
	
	Dim_spectrum=n_elements(lambda_rest)

	print, "        -Interpolate to 1A: done"

	;Create mask ( 10 -> NAN pixels, 6 -> Slicer region)
	mask=intarr(Dim_spectrum)
	mask(*)=0
	nan_index = WHERENAN(spectrum,n_nan) ;Mask NAN regions
		IF n_nan GT 1 then begin 
		spectrum(nan_index)=0.
		error(nan_index)=0.
		mask(nan_index)=10
		ENDIF

	IF KEYWORD_SET(MASK_SLICER) then begin
	slicer=WHERE(lambda_rest GT 5500.*(1+redshift) and lambda_rest LT 6000.*(1+redshift),cnt_slicer ) ;Mask slicer region
	IF cnt_slicer GT 1 then mask(slicer)=6
	ENDIF
	
	IF KEYWORD_SET(MASK_BLUE) then begin
	bleu=WHERE(lambda_rest LT 5500./(1+redshift),cnt_bleu ) ;Mask slicer region
	IF cnt_bleu GT 1 then mask(bleu)=6
	ENDIF


	print, "        -Create mask: done"

	IF KEYWORD_SET(MASK_EMI) then begin
	; Mask emission lines 
	Mask_lines=GAMA_LOAD_emiMask( lambda_rest)	
	mask(WHERE(Mask_lines EQ 1.))=1
	ENDIF




	
	new_struct_spec={Lambda:lambda_rest,Spectrum:spectrum,Error:error,Sensibility:Sensibility,Sky:sky,Mask:mask,name:struct_spec.name, header:struct_spec.header}


return, new_struct_spec
END
