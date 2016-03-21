

;-----------------------------------------------------------------------
; PROCEDURE GAMA_FIT_CONTINUM_LOW_SN
; DESCRIPTION : Fit the stellar continuum with wavelet 
;               for low S_N spectrum
;---------------------------------------------------------------------

FUNCTION GAMA_FIT_CONTINUM_LOW_SN,lambda_rest,spectrum

@GAMA_Main_Config
@GAMA_Config_Sensibility
	
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
		
		continuum=GAMA_WAVELET_CONTINUUM(lambda_rest,spectrum,all_mask,11,5,SILENT=SILENT)
		;oplotObj1 = Obj_New('cgOverPlot',lambda_rest,continuum,thick=2, color="red") 
		;oplotObj2 = Obj_New('cgOverPlot',lambda_rest,all_mask,thick=1, color="green") 
	;	oplotObj3 = Obj_New('cgOverPlot',lambda_rest,Spikes,thick=1, color="pink") 
		
		;cgplot,lambda_rest,smooth(spectrum,3), OPLOTS=[oplotObj1,oplotObj2],xrange=[3700.,6700.]
	
		 return,continuum
END



;-----------------------------------------------------------------------
; PROCEDURE MARVIN_add_config
; DESCRIPTION : Add the parameter from the configuration file in the 
;               starlight grid
;---------------------------------------------------------------------

PRO GAMA_add_config,Comd_config
@GAMA_Config_Sensibility

S_Lambda_MIN_SN=STRTRIM(Lambda_MIN_SN,2)
S_Lambda_MAX_SN=STRTRIM(Lambda_MAX_SN,2)
Comd_pathbase="'s%$PATH_BASE%"+PATH_BASE+"%' "
Comd_pathobs="'s%$PATH_OBS%"+PATH_FIT+"%' "
Comd_pathmask="'s%$PATH_MASK%"+PATH_MASK+"%' "
Comd_pathfit="'s%$PATH_FIT%"+PATH_FIT+"%' "
Comd_red="'s/$REDDENING_LAW/"+REDDENING_LAW+"/' "
Comd_mask="'s/$MASK_STARLIGHT/"+MASK_STARLIGHT+"/' "
Comd_base="'s/$BASE_STARLIGHT/"+BASE_STARLIGHT+"/' "
Comd_configfile="'s/$CONFIG_STARLIGHT/"+CONFIG_STARLIGHT+"/' "
Comd_sn_region1="'s/$Lambda_MIN_SN/"+S_Lambda_MIN_SN+"/' "
Comd_sn_region2="'s/$Lambda_MAX_SN/"+S_Lambda_MAX_SN+"/' "

Comd_config=Comd_red+"-e "+Comd_mask+"-e "+Comd_base+"-e "+Comd_configfile+"-e "+Comd_sn_region1+"-e "+Comd_sn_region2+"-e "+Comd_pathbase+"-e "+Comd_pathobs+"-e "+Comd_pathmask+"-e "+Comd_pathfit
END


;-----------------------------------------------------------------------
; FUNCTION GAMA_SubContinum
; DESCRIPTION : Estimate the stellar continuum in a galaxy spectrum
;               Use the method from Moustakas et al. 2011
; 				3-steps procedure : 
;				  1) Guess v0 and Vd in the bleu region of the spectrum
;				  2) Correct the spectrum from instrumental residuals 
;				  3) Final fit
; INPUT : Spectrum Struct
; RETURN : Error Flag
; SAVE : FITS file raw GAMA spec file + with stellar continuum vector
;---------------------------------------------------------------------


FUNCTION GAMA_Sensibility, struct_spec, PLOT_FIT=PLOT_FIT

@GAMA_Config_Sensibility

spectrum_orig=struct_spec.flux
error_orig=struct_spec.error
sky_orig=struct_spec.sky
lambda=struct_spec.lambda
objet=struct_spec.name


; ---------------------------------------------------
; Step 0 : Preparing the spectrum for Starlight fitting
;       
; ---------------------------------------------------

print, "   [STEP 0] : Preparing spectrum"

	;TRIM if too noisy spectrum at the beginning of the spectrum 
	;relative_error=abs(error_orig/spectrum_orig)
	;mean_error_bleu=mean(relative_error(WHERE( lambda LT 4000)),/NAN)
	;IF mean_error_bleu GT 3. then begin
	;max_value_lamb=max(lambda)- 50.
	;spectrum_orig=GAMA_TRIM_SPECTRUM(lambda,spectrum_orig,[4000.,max_value_lamb])
	;error_orig=GAMA_TRIM_SPECTRUM(lambda,error_orig,[4000.,max_value_lamb])
	;sky_orig=GAMA_TRIM_SPECTRUM(lambda,sky_orig,[4000.,max_value_lamb])
	;lambda=GAMA_TRIM_SPECTRUM(lambda,lambda,[4000.,max_value_lamb])
;	print, "        -Trim spectrum under 4000ang: done"
;	endif	

;Correct from galatic extinction
	;Read value of the extinction from schlegel extinction map	and correct from foreground dust
	;Extinction_filter=[5.155,3.793,2.751,2.086,1.479] ; From schlegel conversion
	RA=SXPAR(struct_spec.header,'RA')
	DEC=SXPAR(struct_spec.header,'DEC')
	euler, RA, DEC, l, b, 1 ;Convert to galatic coordinates
	EBV = dust_getval(l, b)
	CCM_UNRED, lambda, spectrum_orig, EBV
	sxaddpar,struct_spec.header,'EBV',EBV
	
	redshift=SXPAR(struct_spec.header,'Z')
	NAXIS1=SXPAR(struct_spec.header,'NAXIS1')
	lambda_o_rest=lambda/(1+redshift)
	spectrum_orig=spectrum_orig*(1+redshift)
	

	timer_start= SYSTIME(1)  

	;Normalize to Starlight sampling (1A)
	lambda_min=floor(min(lambda_o_rest, /NAN))
	lambda_max=ceil(max(lambda_o_rest,/NAN))
	lambda_rest=indgen(lambda_max-lambda_min)+lambda_min
	spectrum=INTERPOL_V8(spectrum_orig,lambda_o_rest,lambda_rest,/NAN)
	error=INTERPOL_V8(error_orig,lambda_o_rest,lambda_rest,/NAN)
	sky=INTERPOL_V8(sky_orig,lambda_o_rest,lambda_rest,/NAN)
		
	Dim_spectrum=n_elements(spectrum)
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
	slicer=WHERE(lambda_rest GT 5500.*(1+redshift) and lambda_rest LT 6000.*(1+redshift),cnt_slicer ) ;Mask slicer region
	IF cnt_slicer GT 1 then mask(slicer)=6
	print, "        -Create mask: done"
	
	; Mask emission lines 
	Mask_lines=GAMA_LOAD_emiMask( lambda_rest)
	
	; cut in the signal to noise
	SN_from_spectrum=GAMA_SN_continum( lambda_rest, spectrum, [Lambda_MIN_SN,Lambda_MAX_SN])
	IF SN_from_spectrum LT SN_STARLIGHT_CUT then begin
		print, "        !! SN="+TRIM(SN_from_spectrum,2)+" lower than the SN cut (<"+STRCOMPRESS(string(SN_STARLIGHT_CUT),/REMOVE_ALL)+"). Fit Aborted"
		
		return, -4
		ENDIF


; ---------------------------------------------------
; Step 1 : First RUN of Starlight in the blue region
;          to determine vo and vd
; ---------------------------------------------------

print, "   [STEP 1] : Determine vo and vd"

	; Prepare input file for starlight
	file_cxt=objet+sufixe_step1+".cxt"
	openw, U,PATH_FIT+file_cxt ,/GET_LUN
	 	for j=0,Dim_spectrum-1 do begin
			 printf,U,format='(f,e,e,i)',lambda_rest(j),spectrum(j),error(j),mask(j)
		endfor
	FREE_LUN, U 	
	print, "        -Input file: "+file_cxt
	
	; Configure the grid for starlight
	lamnda_range=STRTRIM([lambda_min,lambda_max],2)
	file_out_step1=objet+sufixe_step1+".BN"
	file_grid=objet+sufixe_step1+".in"
	GAMA_add_config,Comd_config
	Comd_lambda1="'s/$Lambda_MIN_FIT/"+lamnda_range(0)+"/' "
	Comd_lambda2="'s/$Lambda_MAX_FIT/4800/' "
	Comd_fitkine="'s/$KIN_FIT/FIT/' "
	Comd_input_file="'s/$INPUT/"+file_cxt+"/' "
	Comd_output_file="'s/$OUTPUT/"+file_out_step1+"/' "
	Comd_vo="'s/$V0_STARLIGHT/0/' "
	Comd_vd="'s/$VD_STARLIGHT/150/' "
	
	command_line="sed -e "+Comd_config+"-e "+Comd_lambda1+"-e "+Comd_lambda2+"-e "+Comd_fitkine+"-e "+Comd_input_file+"-e "+Comd_output_file+"-e "+Comd_vo+"-e "+Comd_vd+ "< grid_file.in > "+ PATH_FIT+file_grid
	spawn,command_line
	print, "        -Input grid: "+file_grid
	
	; Run step 1 in starlight	
	file_log_step1=objet+sufixe_step1+".log"
	start_starlight_step1=PATH_STARLIGHT+"StarlightChains_v04.exe < "+PATH_FIT+file_grid	
	tofile= PATH_FIT+file_log_step1
	spawn,start_starlight_step1,listen
	
	; Read results
	result_fit_exist=file_exist(PATH_FIT+file_out_step1)
		case result_fit_exist of 
		'1': begin
		RDFILE,PATH_FIT+file_out_step1,50,50,line50; chi2
		RDFILE,PATH_FIT+file_out_step1,58,58,line58; v0_min
		RDFILE,PATH_FIT+file_out_step1,59,59,line59; vd_min
		RDFILE,PATH_FIT+file_out_step1,31,31,line31; SN in window
		Chi2_starlight=STRMID(line50,0,20)
		SN_starlight=STRMID(line31,0,20)
		v0_starlight=STRCOMPRESS(STRMID(line58,0,14),/REMOVE_ALL)
		vd_starlight=STRCOMPRESS(STRMID(line59,0,14),/REMOVE_ALL)
		print,"        Starlight fitting: done"
		print, "        -Log file: "+file_log_step1
		print, "        -vo and vd: "+STRCOMPRESS(v0_starlight,/REMOVE_ALL)+" "+STRCOMPRESS(vd_starlight,/REMOVE_ALL)
		print, "        -Chi 2: "+STRCOMPRESS(Chi2_starlight,/REMOVE_ALL)
		;print,"!! SN=",SN_starlight
		end
		'0': begin
		 print, "        !! Fit failed in step 1"
		 return, -1
		end
	endcase	
	
		IF SN_starlight LE SN_STARLIGHT_CUT then begin
		print, "        !! SN"+TRIM(SN_starlight,2)+" lower than the SN cut (<"+STRCOMPRESS(string(SN_STARLIGHT_CUT),/REMOVE_ALL)+"). "
					
		return, -4

		
		ENDIF
	
	


; ---------------------------------------------------
; Step 2 : First RUN of Starlight in the full spectrum
;          to correct from the residuals
; ---------------------------------------------------

	print, "   [STEP 2] : Correct residuals"
	print, "        -Input file: "+file_cxt
	
	;Fit the full spectrum with fix kinematics and mask only emission lines
	file_out=objet+sufixe_step2+".BN"
	file_grid=objet+sufixe_step2+".in"
	GAMA_add_config,Comd_config
	Comd_lambda1="'s/$Lambda_MIN_FIT/"+lamnda_range(0)+"/' "
	Comd_lambda2="'s/$Lambda_MAX_FIT/"+lamnda_range(1)+"/' "
	Comd_fitkine="'s/$KIN_FIT/FXK/' "
	Comd_input_file="'s/$INPUT/"+file_cxt+"/' "
	Comd_output_file="'s/$OUTPUT/"+file_out+"/' "
	Comd_vo="'s/$V0_STARLIGHT/"+v0_starlight+"/' "
	Comd_vd="'s/$VD_STARLIGHT/"+vd_starlight+"/' "	
	command_line="sed -e "+Comd_config+"-e "+Comd_lambda1+"-e "+Comd_lambda2+"-e "+Comd_fitkine+"-e "+Comd_input_file+"-e "+Comd_output_file+"-e "+Comd_vo+"-e "+Comd_vd+ "< grid_file.in > "+ PATH_FIT+file_grid
	spawn,command_line,listen
	print, "        -Input grid: "+file_grid
	
	;Run step 2 fit
	file_log_step2=objet+sufixe_step2+".log"
	start_starlight_step2=PATH_STARLIGHT+"StarlightChains_v04.exe < "+PATH_FIT+file_grid
	spawn,start_starlight_step2,tofile
	
	;Extract results
	result_fit_exist=file_exist(PATH_FIT+file_out)
		case result_fit_exist of 
			'1': begin
			  ;Lecture spectre et FIT from Starlight
			  toto = READ_ASCII(PATH_FIT+file_out,COUNT=n_rows)
			  RDFILE,PATH_FIT+file_out,217,n_rows+17,spectra_out_step2
			  lambda_step2=float(STRMID(spectra_out_step2,1,8))
			  flux_step2=float(STRMID(spectra_out_step2,12,19))
			  fit_step2=float(STRMID(spectra_out_step2,22,29))
			  weight_step2=float(STRMID(spectra_out_step2,32,39))
			  RDFILE,PATH_FIT+file_out,26,26,line26; normalisation
			  normalisation_starlight=float(STRMID(line26,0,16))
	;cgplot,lambda_step2,flux_step2
	;cgoplot,lambda_step2,fit_step2,color='red'
	;pause
			  residual=flux_step2/fit_step2	
			 ; IF KEYWORD_SET(WAVELET) then 
			  sensibility_step2=GAMA_WAVELET_CONTINUUM(lambda_step2,residual,Mask_lines,Wavelet_elements,Sensibilibty_elements,/SILENT)
			  ;sliding_spec=GAMA_SLIDING_MOMENTS(residual,300)	
			 ; sensibility_step2 =sliding_spec[*,0]
		lambda_derest=lambda_rest*(1+redshift)
		lambda_demin=floor(min(lambda_derest, /NAN))
		lambda_demax=ceil(max(lambda_derest,/NAN))
		Dim_spectrum=lambda_demax - lambda_demin
		lambda_new=indgen(Dim_spectrum)+lambda_demin
		print,minmax(lambda_new)
		spectrum_de_z= INTERPOL_V8(spectrum,lambda_derest,lambda_new,/NAN)
		spectrum=INTERPOL(spectrum,lambda_derest,lambda_new)
		error=INTERPOL(error,lambda_derest,lambda_new)
		sensibility=INTERPOL(sensibility_step2,lambda_step2*(1+redshift),lambda_new)
	    ;cgplot,lambda_derest,residual
	    ;cgoplot,lambda_new,sensibility, color='red'

			  print, "        Calculate residuals: done"
			  print, "        Calculate sensibility function : done"
			  
			 
			 
		; Re-add THE NAN value in the spectrum_correct and error
			bad_pixels=WHERE(mask EQ 10,n_bp)
			IF n_bp GE 1 then begin
			spectrum(bad_pixels)=!Values.F_NAN
			error(bad_pixels)=!Values.F_NAN
			endif 
			
			mask(WHERE(Mask_lines EQ 1.))=1
			DATA_final=fltarr(Dim_spectrum,3)
			DATA_final[*,0]=spectrum
			DATA_final[*,1]=error
			DATA_final[*,2]=sensibility

 			;Create spectrum with continous lambda array
			 new_CRVAL1=lambda_demin
			 new_CDELT=1.0
			 new_NAXIS1=Dim_spectrum
			 

			  ;Write in fits file
			header=struct_spec.header
			sxaddpar,header,'CRVAL1',new_CRVAL1
			sxaddpar,header,'NAXIS1',new_NAXIS1
			sxaddpar,header,'CD1_1',new_CDELT
			sxaddpar,header,'EBV',EBV
			sxaddpar,header,'CRPIX1','0'
		    sxaddpar,header,'ROW1',"Spectrum"
		    sxaddpar,header,'ROW2',"Error"
		    sxaddpar,header,'ROW3',"Sensibility"
		    SXDELPAR,header,'ROW4'
		    SXDELPAR,header,'ROW5'
		 	mwrfits,DATA_final,PATH_FIT+preffixe_final+objet+'.fits',header,/Silent,/Create

 			return,1

			end
			'0': begin
			  print, "        !! Fit failed in step 2"
			  return, -2
			  end
	endcase	
	


END


PRO GAMA_Sensibility_ALL

@GAMA_Main_Config
@GAMA_Config_Sensibility

readcol,path_catalogue+catalogue, format="a,l,f,a",spectrumid, cataid, z, QC_Flag
n_obj=n_elements(spectrumid)
QF=intarr(n_obj)
QF(*)=-9 ; not fitted yet

TO_FIT=WHERE( STRCMP(QC_Flag, 'G') EQ 1, count_obj) 

print, "--------------------------------------------- "
print, "-    GAMA Subtract continuum routine  -"
print, "  "
print, " "
print, " Data DIR : "+PATH_OBS
print, "  Working DIR : "+PATH_FIT
print, "  Catalogue loaded :"+catalogue
print, "  Number of spectrum: "+TRIM(n_obj,2)
print, "  Number of spectrum to fit: "+TRIM(count_obj,2)
print, "  Maximum time of computing : "+ TRIM((count_obj*2)/60.,2)+" hours"
print, "--------------------------------------------- "
print, " "
print, " "

;for i=16,16 do begin
for i=0d,count_obj -1 do begin
		file=spectrumid(TO_FIT(i))+".fit"
		print, " "
		print, "[ File "+trim(i+1,2)+"/"+trim(count_obj,2)+"]:"+file
		;print,PATH_OBS+file
		Spec_struc=GAMA_READ_FIT(PATH_OBS+file)
		check_data=size(Spec_struc)
				IF (check_data(0) EQ 0) then begin 
				print, "        !! Failed to open file "+file
				QF(TO_FIT(i))=0	
				ENDIF ELSE begin 
				cp_comd="cp "+PATH_OBS+file+" "+PATH_FIT+file
				QF(TO_FIT(i))=GAMA_Sensibility(Spec_struc)
				ENDELSE	
				
		IF ((i MOD AUTO_SAVE_P) EQ 0) then begin
		 save, spectrumid,cataid,z, QC_Flag,QF, filename=PATH_FIT+saved_file
		ENDIF
	
	endfor
print, "Feed the baby !"
END
