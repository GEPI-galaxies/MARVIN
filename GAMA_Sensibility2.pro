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



FUNCTION GAMA_Sensibility2, struct_spec, PLOT_FIT=PLOT_FIT

@GAMA_Main_Config

objet=struct_spec.name




	spectrum_str=GAMA_PREPARE_STARLIGHT( struct_spec, /EXTINCTION, /MASK_SLICER,/MASK_BLUE,/TRIM)
	objet=spectrum_str.name
	Dim_spectrum=n_elements(spectrum_str.Lambda)	
	
	SN_from_spectrum=GAMA_SN_continum(spectrum_str.Lambda, spectrum_str.spectrum, [Lambda_MIN_SN,Lambda_MAX_SN])
	IF SN_from_spectrum LT SN_STARLIGHT_CUT then begin
		print, "        !! SN="+TRIM(SN_from_spectrum,2)+" lower than the SN cut (<"+STRCOMPRESS(string(SN_STARLIGHT_CUT),/REMOVE_ALL)+"). Fit Aborted"
		
		return, -4
	ENDIF


	fit_flag=GAMA_RUN_STARLIGHT(objet, spectrum_str,sufixe_step1)
	
	
	; Read results
	case fit_flag of 
			'1': begin
			  file_out=objet+sufixe_step1+".BN"
			  ;Lecture spectre et FIT from Starlight
			  toto = READ_ASCII(PATH_FIT+file_out,COUNT=n_rows)
			  RDFILE,PATH_FIT+file_out,217,n_rows+17,spectra_out_step2
			  lambda_step1=float(STRMID(spectra_out_step2,1,8))
			  flux_step1=float(STRMID(spectra_out_step2,12,19))
			  fit_step1=float(STRMID(spectra_out_step2,22,29))
			  weight_step1=float(STRMID(spectra_out_step2,32,39))
			  RDFILE,PATH_FIT+file_out,26,26,line26; normalisation
			  normalisation_starlight=float(STRMID(line26,0,16))

			  Residual=flux_step1/fit_step1	
			  Mask_lines=GAMA_LOAD_emiMask(lambda_step1)	
			  ;Sensibility=GAMA_WAVELET_CONTINUUM(lambda_step1,Residual,Mask_lines,Wavelet_elements_SENS,Sensibilibty_elements_SENS,/SILENT)
			Sensibility= GAMA_SMOOTH_CONTINUUM(lambda_step1,Residual,Mask_lines,11)

			  print, "        Calculate sensibility correction: done"
			  Sensibility=INTERPOL(Sensibility,lambda_step1,spectrum_str.Lambda)

			  print, "        Correct spectrum: done"
			  end
			'0': begin
			  print, "        !! Fit failed Sensibility correction"
			  return, -2
			  end
		endcase	

		; Re-add THE NAN value in the spectrum_correct and error
		bad_pixels=WHERE(spectrum_str.Mask EQ 10,n_bp)
		IF n_bp GE 1 then begin
		spectrum_str.Spectrum(bad_pixels)=!Values.F_NAN
		spectrum_str.error(bad_pixels)=!Values.F_NAN
		ENDIF 
		
		; Add emission lines to mask 
		Mask_lines=GAMA_LOAD_emiMask(spectrum_str.Lambda)	
		spectrum_str.Mask(WHERE(Mask_lines EQ 1.))=1
		
		;De-redshift and remove remove the dust correction
		EBV= SXPAR(struct_spec.header,'EBV')
		redshift=SXPAR(struct_spec.header,'Z')
		CCM_UNRED, spectrum_str.Lambda, spectrum_str.Spectrum , (-1.)*EBV
		spectrum_str.Lambda=spectrum_str.Lambda*(1+redshift)
		spectrum_str.Spectrum= spectrum_str.Spectrum/(1+redshift)

		;Create spectrum with continous lambda array
		 new_CRVAL1=spectrum_str.Lambda(0)
		 new_CDELT=1.0*(1+redshift)
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
	  	sxaddpar,header,'ROW4',"Sky"
	 	sxaddpar,header,'ROW5',"Mask"

		lambda_temp=indgen(new_NAXIS1)*new_CDELT+new_CRVAL1
	
	
		DATA_final=fltarr(Dim_spectrum,6)
		DATA_final[*,0]=spectrum_str.Spectrum
		DATA_final[*,1]=spectrum_str.Error
		DATA_final[*,2]=Sensibility
		DATA_final[*,3]=spectrum_str.Sky
		DATA_final[*,4]=spectrum_str.Mask
		WRITEFITS,PATH_SENS+'Fit/Sens_'+objet+'.fits',DATA_final,header
 		cgplot, lambda_step1*(1+redshift), smooth(Residual,6),yrange=[0,2],xrange=[3700,8900]
		cgoplot,lambda_temp,Sensibility,color='green'

	END
	


PRO GAMA_Sensibility_ALL2

@GAMA_Main_Config

readcol,PATH_CATALOGUE+CATALOGUE_SENS_INPUT, format="a,l,f",spectrumid, cataid, z
n_obj=n_elements(spectrumid)
QF=intarr(n_obj)
QF(*)=-9 ; not fitted yet

print, "--------------------------------------------- "
print, "-    GAMA Subtract continuum routine  -"
print, "  "
print, " "
print, " Data DIR : "+PATH_SENS+'Data/'
print, "  Working DIR : "+PATH_SENS+'Fit/'
print, "  Catalogue loaded :"+CATALOGUE_SENS_INPUT
print, "  Number of spectrum: "+TRIM(n_obj,2)
print, "  Maximum time of computing : "+ TRIM((n_obj*1.)/60.,2)+" hours"
print, "--------------------------------------------- "
print, " "
print, " "

;for i=16,16 do begin
for i=0d, n_obj -1 do begin
		file=spectrumid(i)+".fit"
		print, " "
		print, "[ File "+trim(i+1,2)+"/"+trim(count_obj,2)+"]:"+file
		Spec_struc=GAMA_READ_FITS(PATH_SENS+"/Data/"+file)
		check_data=size(Spec_struc)
				IF (check_data(0) EQ 0) then begin 
				print, "        !! Failed to open file "+file
				QF(i)=0	
				ENDIF ELSE begin 
				QF(i)=GAMA_Sensibility2(Spec_struc)
				ENDELSE	
				
		IF ((i MOD AUTO_SAVE_SENS) EQ 0) then begin
		 save, spectrumid,cataid,z,QF, filename=PATH_SENS+'Sens_'+CATALOGUE_SENS_INPUT+'.sav'
		ENDIF
	
	endfor
print, "Feed the baby !"
END	

