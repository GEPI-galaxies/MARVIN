
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


FUNCTION GAMA_FitContinum, struct_spec, PLOT_FIT=PLOT_FIT

@GAMA_Main_Config

; ---------------------------------------------------
; Step 0 : Prepare the spectrum for STARLIGHT
;          Test the spectrum S/N. If bellow the SNR threshold, 
;	   estimate the continum using a smoothing 
; ---------------------------------------------------

	spectrum_str=GAMA_PREPARE_STARLIGHT( struct_spec, /MASK_SLICER) ;, /TRIM, /EXTINCTION
	objet=spectrum_str.name
	Dim_spectrum=n_elements(spectrum_str.Lambda)

	; Verify the signal to noise in the continum. If bellow the thresold evaluate the continum with smoothing
	SN_from_spectrum=GAMA_SN_continum( spectrum_str.Lambda, spectrum_str.Spectrum, [Lambda_MIN_SN,Lambda_MAX_SN])
		
		IF SN_from_spectrum LT SN_STARLIGHT_CUT then begin
		print, "        !! SN="+TRIM(SN_from_spectrum,2)+" lower than the SN cut (<" $
					+STRCOMPRESS(string(SN_STARLIGHT_CUT),/REMOVE_ALL)+"). Fit Aborted"
		continuum = GAMA_FIT_CONTINUM_LOW_SN(spectrum_str.Lambda, spectrum_str.Spectrum)
			
			IF ~TAG_EXIST(spectrum_str,'Sensibility') then begin
			Sensibility=fltarr(Dim_spectrum)
			Sensibility(*)=1.
			ENDIF

			IF TAG_EXIST(spectrum_str,'Sensibility') then begin
			Sensibility=spectrum_str.Sensibility
			ENDIF



		;Create spectrum with continous lambda array
		 new_CRVAL1=spectrum_str.Lambda(0)
		 new_CDELT=1.0
		 new_NAXIS1=Dim_spectrum
		 
		;Update the new header
		header=spectrum_str.header
		sxaddpar,header,'CRVAL1',new_CRVAL1
		sxaddpar,header,'NAXIS1',new_NAXIS1
		sxaddpar,header,'CD1_1',new_CDELT
		sxaddpar,header,'CRPIX1','0'
	  	sxaddpar,header,'ROW1',"Spectrum"
	   	 sxaddpar,header,'ROW2',"Error"
	    	sxaddpar,header,'ROW3',"Sensibility"
	   	 sxaddpar,header,'ROW4',"Continuum"
	 	sxaddpar,header,'ROW5',"Sky"
	 	sxaddpar,header,'ROW6',"Mask"
		sxaddpar,header,'FITC',"Wavelet"

		; Re-add THE NAN value in the spectrum_correct and error
		bad_pixels=WHERE(spectrum_str.Mask EQ 10,n_bp)
		IF n_bp GE 1 then begin
		spectrum_str.Spectrum(bad_pixels)=!Values.F_NAN
		spectrum_str.Error(bad_pixels)=!Values.F_NAN
		endif 
		
		; Add emission lines to mask 
		Mask_lines=GAMA_LOAD_emiMask(spectrum_str.Lambda)	
		spectrum_str.Mask(WHERE(Mask_lines EQ 1.))=1
		

		DATA_final=fltarr(Dim_spectrum,6)
		DATA_final[*,0]=spectrum_str.Spectrum
		DATA_final[*,1]=spectrum_str.error
		DATA_final[*,2]=Sensibility
		DATA_final[*,3]=continuum
		DATA_final[*,4]=spectrum_str.Sky
		DATA_final[*,5]=spectrum_str.Mask

		IF KEYWORD_SET(PLOT_FIT) then begin
		oplotObj1 = Obj_New('cgOverPlot',spectrum_str.Lambda ,continuum,thick=1, color="red") 
		oplotObj2 = Obj_New('cgOverPlot',spectrum_str.Lambda,smooth(spectrum_str.Spectrum,6),thick=1) 
		cgplot,spectrum_str.Lambda ,smooth(spectrum_str.Spectrum*Sensibility,6), OPLOTS=[oplotObj1,oplotObj2], $
						 color="grey",ytitle='Flux [10^-16 erg/s/cm^2/A] ',xtitle='Wavelength [A]',/YLOG,yrange=[0.1,max(smooth(spectrum_str.Spectrum,6))]
		ENDIF	
		WRITEFITS,PATH_FIT+preffixe_final+objet+'.fits',DATA_final,header

		print, "        -Evaluate continum with wavelet: done"

		return, -4
		ENDIF



; ---------------------------------------------------
; Step 1 : First RUN of Starlight in the blue region
;          to determine vo and vd
; ---------------------------------------------------

print, "   [STEP 1] : Determine vo and vd"
	lambda_min=min(spectrum_str.Lambda)
	fit_flag=GAMA_RUN_STARLIGHT(objet, spectrum_str,lambda_vector=[lambda_min,4800],sufixe_step1)
	
	; Read results
	case fit_flag of 
		'1': begin
		RDFILE,PATH_FIT+objet+sufixe_step1+".BN",50,50,line50; chi2
		RDFILE,PATH_FIT+objet+sufixe_step1+".BN",58,58,line58; v0_min
		RDFILE,PATH_FIT+objet+sufixe_step1+".BN",59,59,line59; vd_min
		RDFILE,PATH_FIT+objet+sufixe_step1+".BN",31,31,line31; SN in window
		Chi2_starlight=STRMID(line50,0,20)
		SN_starlight=STRMID(line31,0,20)
		v0_starlight=STRCOMPRESS(STRMID(line58,0,14),/REMOVE_ALL)
		vd_starlight=STRCOMPRESS(STRMID(line59,0,14),/REMOVE_ALL)
		print, "        -Log file: "+ objet+sufixe_step1+".log"
		print, "        -vo and vd: "+STRCOMPRESS(v0_starlight,/REMOVE_ALL)+" "+STRCOMPRESS(vd_starlight,/REMOVE_ALL)
		print, "        -Chi 2: "+STRCOMPRESS(Chi2_starlight,/REMOVE_ALL)
		end
		'0': begin
		 print, "        !! Fit failed in step 1"
		 return, -1
		end
	endcase	
	
		IF SN_starlight LE SN_STARLIGHT_CUT then begin
		print, "        !! SN"+TRIM(SN_starlight,2)+" lower than the SN cut (<"+STRCOMPRESS(string(SN_STARLIGHT_CUT),/REMOVE_ALL)+"). "
			
		continuum = GAMA_FIT_CONTINUM_LOW_SN(spectrum_str.Lambda,spectrum_str.Spectrum)	

			IF ~TAG_EXIST(spectrum_str,'Sensibility') then begin
			Sensibility=fltarr(Dim_spectrum)
			Sensibility(*)=1.
			ENDIF

			IF TAG_EXIST(spectrum_str,'Sensibility') then begin
			Sensibility=spectrum_str.Sensibility
			ENDIF
		
		;Create spectrum with continous lambda array
		 new_CRVAL1=spectrum_str.Lambda(0)
		 new_CDELT=1.0
		 new_NAXIS1=Dim_spectrum
		 
		;Update the new header
		header=spectrum_str.header
		sxaddpar,header,'CRVAL1',new_CRVAL1
		sxaddpar,header,'NAXIS1',new_NAXIS1
		sxaddpar,header,'CD1_1',new_CDELT
		sxaddpar,header,'CRPIX1','0'
	  	sxaddpar,header,'ROW1',"Spectrum"
	        sxaddpar,header,'ROW2',"Error"
	        sxaddpar,header,'ROW3',"Sensibility"
	        sxaddpar,header,'ROW4',"Continuum"
	  	sxaddpar,header,'ROW5',"Sky"
	 	sxaddpar,header,'ROW6',"Mask"
		sxaddpar,header,'FITC',"Wavelet"
		
	; Re-add THE NAN value in the spectrum_correct and error
		bad_pixels=WHERE(spectrum_str.Mask EQ 10,n_bp)
		IF n_bp GE 1 then begin
		spectrum_str.Spectrum(bad_pixels)=!Values.F_NAN
		spectrum_str.error(bad_pixels)=!Values.F_NAN
		endif 
		
		; Add emission lines to mask 
		Mask_lines=GAMA_LOAD_emiMask(spectrum_str.Lambda)	
		spectrum_str.Mask(WHERE(Mask_lines EQ 1.))=1
	
		DATA_final=fltarr(Dim_spectrum,6)
		DATA_final[*,0]=spectrum_str.Spectrum
		DATA_final[*,1]=spectrum_str.Error
		DATA_final[*,2]=Sensibility
		DATA_final[*,3]=continuum
		DATA_final[*,4]=spectrum_str.Sky
		DATA_final[*,5]=spectrum_str.Mask
		WRITEFITS,PATH_FIT+preffixe_final+objet+'.fits',DATA_final,header
		
		IF KEYWORD_SET(PLOT_FIT) then begin
		oplotObj1 = Obj_New('cgOverPlot',spectrum_str.Lambda ,continuum,thick=1, color="red") 
		oplotObj2 = Obj_New('cgOverPlot',spectrum_str.Lambda,smooth(spectrum_str.Spectrum,6),thick=1) 
		cgplot,spectrum_str.Lambda ,smooth(spectrum_str.Spectrum*Sensibility,6), OPLOTS=[oplotObj1,oplotObj2], $
						 color="grey",ytitle='Flux [10^-16 erg/s/cm^2/A] ',xtitle='Wavelength [A]',/YLOG,yrange=[0.1,max(smooth(spectrum_str.Spectrum,6))]
		ENDIF	

		print, "        -Evaluate continum with wavelet: done"
		return, -4
		ENDIF
	

; ---------------------------------------------------
; Step 2 : First RUN of Starlight in the full spectrum
;          to correct from the residuals
; ---------------------------------------------------

	print, "   [STEP 2] : Correct residuals"
	Kinematics_vector_step1=[float(STRCOMPRESS(v0_starlight,/REMOVE_ALL)),float(STRCOMPRESS(vd_starlight,/REMOVE_ALL))]
	fit_flag=GAMA_RUN_STARLIGHT(objet, spectrum_str,sufixe_step2,Kinematics_vector=Kinematics_vector_step1,/FIX_KINEMATICS)
	
	; Read results
	case fit_flag of 
			'1': begin
			  file_out=objet+sufixe_step2+".BN"
			  ;Lecture spectre et FIT from Starlight
			  toto = READ_ASCII(PATH_FIT+file_out,COUNT=n_rows)
			  RDFILE,PATH_FIT+file_out,217,n_rows+17,spectra_out_step2
			  lambda_step2=float(STRMID(spectra_out_step2,1,8))
			  flux_step2=float(STRMID(spectra_out_step2,12,19))
			  fit_step2=float(STRMID(spectra_out_step2,22,29))
			  weight_step2=float(STRMID(spectra_out_step2,32,39))
			  RDFILE,PATH_FIT+file_out,26,26,line26; normalisation
			  normalisation_starlight=float(STRMID(line26,0,16))

			  residual=flux_step2-fit_step2	
			  Mask_lines=GAMA_LOAD_emiMask(lambda_step2)	
			  residual_step2=GAMA_WAVELET_CONTINUUM(lambda_step2,residual,Mask_lines,Wavelet_elements,Sensibilibty_elements,/SILENT)

			  print, "        Calculate residuals: done"
			  residual=INTERPOL(residual_step2,lambda_step2,spectrum_str.Lambda)*replicate(normalisation_starlight,n_elements(spectrum_str.Lambda))
			  old_spectrum=spectrum_str.Spectrum
			  old_error=spectrum_str.Error
			  spectrum_str.Spectrum = spectrum_str.Spectrum-residual
			  spectrum_str.error=sqrt((abs(residual))^2+spectrum_str.error^2)
			  flux_calib_error=TOTAL(abs((spectrum_str.Spectrum-old_spectrum)/old_spectrum),/NAN)/n_elements(spectrum_str.Spectrum)
			  plot,spectrum_str.Lambda,residual
			  print, "        Correct spectrum: done"
			  print, "        - Residuals correction: "+TRIM(flux_calib_error*100.,2)+"%"
			  end
			'0': begin
			  print, "        !! Fit failed in step 2"
			  return, -2
			  end
	endcase	


; ---------------------------------------------------
; Step 3 : Final RUN of Starlight in the full spectrum
;         
; ---------------------------------------------------


print, "   [STEP 3] : Final Starlight fit"

	
	; Prepare input file for starlight
	un_mask= (WHERE(spectrum_str.Mask EQ 6, n_masked) )
	IF n_masked GE 1 then  spectrum_str.Mask(un_mask)=0 ;remove slicer region from the mask, and keep NAN pixel in the mask
	
	fit_flag=GAMA_RUN_STARLIGHT(objet, spectrum_str,sufixe_step3,Kinematics_vector=Kinematics_vector_step1)
	; Read results
	case fit_flag of 
			'1': begin
	
	;Lecture spectre et FIT from Starlight
	file_out=objet+sufixe_step3+".BN"
	toto = READ_ASCII( PATH_FIT+file_out,COUNT=n_rows)
	RDFILE,PATH_FIT+file_out,217,n_rows+17,spectra_out_step3
	lambda_step3=float(STRMID(spectra_out_step3,1,8))
	flux_step3=float(STRMID(spectra_out_step3,12,19))
	fit_step3=float(STRMID(spectra_out_step3,22,29))
	weight_step3=float(STRMID(spectra_out_step3,32,39))
	RDFILE,PATH_FIT+file_out,50,50,line50; chi2
	RDFILE,PATH_FIT+file_out,58,58,line58; v0_min
	RDFILE,PATH_FIT+file_out,59,59,line59; vd_min
	RDFILE,PATH_FIT+file_out,26,26,line26; normalisation
	Chi2_starlight=STRMID(line50,0,20) 
	normalisation_starlight=float(STRMID(line26,0,16))
	v0_starlight=STRCOMPRESS(STRMID(line58,0,14),/REMOVE_ALL)
	vd_starlight=STRCOMPRESS(STRMID(line59,0,14),/REMOVE_ALL)
	
	fit_normalized=fit_step3*replicate(normalisation_starlight,Dim_spectrum)
	print,"        Starlight fitting: done"
	print, "        -Log file: "+objet+sufixe_step3+".log"
	print, "        -vo and vd: "+STRCOMPRESS(v0_starlight,/REMOVE_ALL)+" "+STRCOMPRESS(vd_starlight,/REMOVE_ALL)
	print, "        -Chi 2: "+STRCOMPRESS(Chi2_starlight,/REMOVE_ALL)
	
	print, "   [FIT FINISH] "
	
	IF KEYWORD_SET(PLOT_FIT) then begin
	oplotObj1 = Obj_New('cgOverPlot',lambda_step3,fit_normalized,thick=1, color="red") 
	oplotObj2 = Obj_New('cgOverPlot',lambda_step3,smooth(spectrum_str.Spectrum,6),thick=1) 
	oplotObj3 = Obj_New('cgOverPlot',lambda_step3,spectrum_str.Error,thick=1, color="violet") 
	cgplot,spectrum_str.Lambda ,smooth(old_spectrum,6), OPLOTS=[oplotObj1,oplotObj2], $
							 color="grey",ytitle='Flux [10^-16 erg/s/cm^2/A] ',xtitle='Wavelength [A]',/YLOG,yrange=[0.1,max(smooth(old_spectrum,6))]
	ENDIF	
	
	IF ~TAG_EXIST(spectrum_str,'sensibility') then begin
	Sensibility=fltarr(Dim_spectrum)
	Sensibility(*)=1.
	ENDIF

	IF TAG_EXIST(spectrum_str,'sensibility') then begin
	Sensibility=spectrum_str.Sensibility
	ENDIF

	 ;Create spectrum with continous lambda array
	 new_CRVAL1=lambda_min
	 new_CDELT=1.0
	 new_NAXIS1=Dim_spectrum
	 
	;Update the new header
		header=spectrum_str.header
		sxaddpar,header,'CRVAL1',new_CRVAL1
		sxaddpar,header,'NAXIS1',new_NAXIS1
		sxaddpar,header,'CD1_1',new_CDELT
		sxaddpar,header,'CRPIX1','0'
		sxaddpar,header,'CRPIX1','0'
	    sxaddpar,header,'ROW1',"Spectrum"
	    sxaddpar,header,'ROW2',"Error"
	    sxaddpar,header,'ROW3',"Sensibility"
	    sxaddpar,header,'ROW4',"Continuum"
	 	sxaddpar,header,'ROW5',"Sky"
	 	sxaddpar,header,'ROW6',"Mask"
		sxaddpar,header,'ROW7',"Residuals"
		sxaddpar,header,'FITC',"Starlight"

	; Re-add THE NAN value in the spectrum_correct and error
	bad_pixels=WHERE(spectrum_str.Mask EQ 10,n_bp)
		IF n_bp GE 1 then begin
		spectrum_str.Spectrum(bad_pixels)=!Values.F_NAN
		spectrum_str.error(bad_pixels)=!Values.F_NAN
		endif 
		
	spectrum_str.Mask(WHERE(Mask_lines EQ 1.))=1
	DATA_final=fltarr(Dim_spectrum,7)
	DATA_final[*,0]=spectrum_str.Spectrum
	DATA_final[*,1]=spectrum_str.Error
	DATA_final[*,2]=Sensibility
	DATA_final[*,3]=fit_normalized
	DATA_final[*,4]=spectrum_str.Sky
	DATA_final[*,5]=spectrum_str.Mask
	DATA_final[*,6]=residual
	WRITEFITS,PATH_FIT+preffixe_final+objet+'.fits',DATA_final,header
	return,1
	end
	
	'0': begin
		 print, "        !! Fit failed in step 3"
		 return, -3
		end
	endcase	
	


END


PRO GAMA_SubContinum_ALL

@GAMA_Main_Config


readcol,PATH_CATALOGUE+CATALOGUE_INPUT_FIT, format="a,l,f,a",spectrumid, cataid, z, QC_Flag
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
print, "  Catalogue loaded :"+CATALOGUE_INPUT_FIT
print, "  Number of spectrum: "+TRIM(n_obj,2)
print, "  Number of spectrum to fit: "+TRIM(count_obj,2)
print, "  Maximum time of computing : "+ TRIM((count_obj*2)/60.,2)+" hours"
print, "--------------------------------------------- "
print, " "
print, " "


for i=0,count_obj -1 do begin
		file=spectrumid(TO_FIT(i))+".fits"
		print, " "
		print, "[ File "+trim(i+1,2)+"/"+trim(count_obj,2)+"]:"+file
		Spec_struc=GAMA_READ_FITS(PATH_OBS+file)
		check_data=size(Spec_struc)
				IF (check_data(0) EQ 0) then begin 
				print, "        !! Failed to open file "+file
				QF(TO_FIT(i))=0	
				ENDIF ELSE begin 
				;cp_comd="cp "+PATH_OBS+file+" "+PATH_FIT+file
				QF(TO_FIT(i))=GAMA_FitContinum(Spec_struc,/PLOT_FIT)
				ENDELSE	
				
		IF ((i MOD AUTO_SAVE_FIT) EQ 0) then begin
		 save, spectrumid,cataid,z, QC_Flag,QF, filename=SAVED_FIT
		ENDIF
	
	endfor
print, "Feed the baby !"
END
