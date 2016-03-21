
pro LineMeasure__define
 
;  This routine defines the structure for the flux measurement grid
 
  tmp = {LineMeasure, $
         Name : '', $			    ; Name of the line
         Center: -999.d,$           ; Centroid of the line
         ErrCenter: -999.d,   $     ; Error Centroid
         Width: -999.d, $           ; Width line
         ErrWidth: -999.d,   $      ; Err Width
         Peak: -999.d, $            ; Peak of the line
         ErrPeak: -999.d,   $       ; Error Peak
         Background: -999.d,   $    ; Background under the line
         ErrBackground: -999.d,   $ ; Error on the background under the line
         Chi2: -999.d,   $          ; Chi2 of the fit
         QC_index:-99, $            ; Quality Flag
         
         ;Derivate measure
         Flux_gaussian: -999.d, $     ;Flux from model the fitted model
		 Err_Flux_gaussian: -999.d, $ ; Error on the estimated Flux of the gaussian
		 Flux_integrated: -999.d, $     ;Flux integrated under the line
		 Err_Flux_integrated: -999.d, $ ; Error on the integrated Flux 
		 Flux: -999.d, $
		 Flux_limit: -999.d, $ ;Flux from model the fitted model
   		 Noise_level: -999.d, $;Noise level around the line
		 SNR: -999.d, $ ;Signal to noise of detection of the line
		 EW: -999.d $ ;EW the line
		 
         }
 
end


	;-----------------------------------------------------------------------
	; FUNCTION Error_Flux_GAUSS
	; DESCRIPTION : Calculate the error on the flux using Monte-Carlo
	;               from the error on the gaussian parameters (sigma, peak)
	;---------------------------------------------------------------------
	
	FUNCTION Error_Flux_GAUSS, Gauss_P, Err_P
	
	N_MC=1000
	flux_in_line=fltarr(N_MC)
	
	For i=0,N_MC-1 do begin 
	peak_fit=Gauss_P[2]+randomn(seed,1)*Err_P[2]
	sigma_fit=Gauss_P[1]+randomn(seed,1)*Err_P[1]
	flux_in_line[i]=peak_fit*sigma_fit*SQRT(2*!DPI);*factor_normalization
	endfor
	err_flux=sigma(flux_in_line)
	return, err_flux
	END
	
	


	;------------------------------------------------------------------------------------------
	; FUNCTION GAMA_multi_fit  
	;
	; DESCRIPTION : Routine to fit simultaneously multiple lines using gaussian profiles
	; 
	; INPUT :  Lambda -> Wavelenght array in Ang 
	;		   Spectrum -> Flux density array (normalized unit)
	;          Mask -> Bad-pixel map array ( 0 Good, >1 Bad)
	;          lines_ref -> array with the wavelenght of reference of the lines to fit 
	;
	; OPTION : /PLOT_FIT  -> Plot fit of the lines
	;		   /FIX_SIGMA -> Fix the sigma of all the gaussians to the first line of the fitting group
	;		   /FIX_separation -> Fix the separation between lines relatively to the first line of the fitting group 
	;							 according to the reference value of the lines
	;			/FIX_ratio between lines in a doublet 
	; 			continuum -> fit EW are measure give the continuum spectrum 
	; CONFIGURATION : Use a configuration file @GAMA_Config_FitLines
	; 		The configuration file contain Fitting parameters read by this function.
	; 		$default_Line_window=80 -> in [A], define the window from both side of the line(s) to consider when fitting
	;							;Should be large enough to cover both line and background to calculate the noise level
	;       $SIGMA_MIN -> in [A], Minimun possible value of the sigma for the gaussian model
	;       $SIGMA_MAX -> in [A], Maxinum possible value of the sigma for the gaussian model
	; 		$SIGMA_DEFAULT -> in [A], Starting value for the sigma value of the gaussian model
	;		$WAVELENGHT_TOLERANCE -> in [A], Tolerance on the value of the gaussian center relative to the reference value
	;
	; OUTPUT : Structure with N_lines X structure LineMeasure (DIM=17 x N_lines)
	;RETURN_VALUES array with DIM=10 x N_lines
	; 		   RETURN_VALUES[0,*] -> center of the gaussian
	; 		   RETURN_VALUES[1,*] -> Sigma of the gaussian
	; 		   RETURN_VALUES[2,*] -> Flux of the gaussian
	; 		   RETURN_VALUES[3,*] -> Error on the estimated Flux of the gaussian
	; 		   RETURN_VALUES[4,*] -> Noise level around the line
	; 		   RETURN_VALUES[5,*] -> Signal to noise of detection of the line
	; 		   RETURN_VALUES[6,*] -> Quality Flag for the gaussian fit - see FLAG_FIT 
	; 		   RETURN_VALUES[7,*] -> Peak amplitude of the line
	; 		   RETURN_VALUES[8,*] -> Continuum fit P[0]
	; 		   RETURN_VALUES[9,*] -> Continuum fit P[1]
	; 
	; FLAG_FIT : This flag gives an error flag for the fitting of a emission line. The possible values are:
	; 				 1 -> Emission line fitted and detected 
	;				-1 -> The noise level couldn't be evaluated. Fail
	; 				-2 -> Relative error on the flux of the lines is above threshold $MAX_error_flux
	;				-3 -> SN of detection under the detection threshold $SN_detection_threshold
	;		-16 to -18 -> MPFIT error 
	;-------------------------------------------------------------------------------------------

    FUNCTION GAMA_multi_gaussian,lambda_line,flux_line,error,mask,lines_ref, $ 
				PLOT_FIT=PLOT_FIT,FIX_SIGMA=FIX_SIGMA,FIX_separation=FIX_separation, FIX_RATIO=FIX_RATIO, ratio=ratio, continuum_line=continuum_line
	
    @GAMA_Main_Config
		
      N_gauss=n_elements(lines_ref)
      centers=fltarr(N_gauss)
      Peak_line=fltarr(N_gauss)
	
      Lines_measurements= REPLICATE({LineMeasure}, N_gauss)
	
      ;Initialize the measure structure    
      Lines_measurements[*].Name=''
      Lines_measurements[*].Center=-999.
      Lines_measurements[*].ErrCenter=-999.
      Lines_measurements[*].Width=-999.
      Lines_measurements[*].ErrWidth=-999.
      Lines_measurements[*].Peak=-999.
      Lines_measurements[*].ErrPeak=-999.
      Lines_measurements[*].Background=-999.
      Lines_measurements[*].ErrBackground=-999.
      Lines_measurements[*].Chi2=-999.
      Lines_measurements[*].QC_index=-999.
      Lines_measurements[*].Flux_gaussian=-999.
      Lines_measurements[*].Err_Flux_gaussian=-999.
      Lines_measurements[*].Flux_limit=-999.
      Lines_measurements[*].Flux_integrated=-999.
      Lines_measurements[*].Err_Flux_integrated=-999.
      Lines_measurements[*].Noise_level=-999.
      Lines_measurements[*].SNR=-999.
      Lines_measurements[*].Flux=-999.
      Lines_measurements[*].EW=-999.
      	 
      ;Define fitting function
        expr='P[0]'
	parinfo = replicate({value:0.D, tied:" ", fixed:0, limited:[0,0], $
	                       limits:[0.D,0]}, N_gauss*3+1)
	start=[0.]
	parinfo[0].limited[0] = 1
	parinfo[0].limits[0]  = -1

	  FOR i=1, N_gauss do begin
		expr=expr+'+ GAUSS1(X, P['+STRTRIM(i*3-2,1)+":"+STRTRIM(i*3,1)+'],/PEAK)'
		index_center=3*i-2
		index_sigma=3*i-1
		index_peak=3*i
		
		;Start Values
		centers[i-1]=lines_ref[i-1]
		Peak_line[i-1]=MAX(flux_line(WHERE(lambda_line GE centers[i-1]-5. and lambda_line LT centers[i-1] +5.)),/NAN)
		start=[start,centers[i-1],SIGMA_DEFAULT,Peak_line[i-1]]

		IF KEYWORD_SET(FIX_separation) and i EQ 2  THEN BEGIN 
	    	parinfo[index_center].tied='p[1]+('+string(lines_ref[0]-lines_ref[i-1])+")"
	        ENDIF ELSE begin
	    	parinfo[index_center].limited[0] = 1
			parinfo[index_center].limits[0]  =  centers[i-1]-WAVELENGTH_TOLERANCE
			parinfo[index_center].limited[1] = 1
			parinfo[index_center].limits[1]  = centers[i-1]+WAVELENGTH_TOLERANCE
		ENDELSE
				
		IF KEYWORD_SET(FIX_RATIO) and i EQ 2 THEN BEGIN 
		parinfo[index_peak].tied  = 'p[3]/ratio'
	        ENDIF ELSE begin
		parinfo[index_peak].limited[0] = 1
		parinfo[index_peak].limits[0]  = 0.D
		ENDELSE
		
		parinfo[index_sigma].limited[0] = 1
		parinfo[index_sigma].limits[0]  = SIGMA_MIN
		parinfo[index_sigma].limited[1] = 1
		parinfo[index_sigma].limits[1]  = 3* SIGMA_MAX
		

	   
	        IF KEYWORD_SET(FIX_SIGMA) THEN BEGIN 
	  		IF i GT 1 then begin
	  		parinfo[index_sigma].tied = 'p[2]'
	  		ENDIF
	        ENDIF 
	ENDFOR
	  
    ;Fit the lines
	parinfo[*].value = start
 	P = mpfitexpr(expr, lambda_line, flux_line, error,yfit=gauss_line,PARINFO=parinfo,$
	           BESTNORM=BESTNORM, PERROR=PERROR,/NAN,STATUS=ERROR_STATUS,/QUIET,DOF=DOF,$
	             best_fjac=best_fjac, pfree_index=pfree_index, /calc_fjac, covar=pcovar)

	IF ERROR_STATUS LT 1 then begin ;Catching error during fit
	Lines_measurements.QC_index = ERROR_STATUS - 20.
	return,Lines_measurements
	endif 
	
	;Estimates the error of each parameter and the error on the flux
	 PCERROR = PERROR * SQRT(BESTNORM / DOF) ; scaled uncertainties
	 ;YCOVAR = mpproperr(best_fjac, pcovar, pfree_index=pfree_index)
	 ;PCERROR = sqrt(YCOVAR) 

    ;Measure properties of the lines	
	  FOR i=1,N_gauss do begin
		Lines_measurements[i-1].Name=''
		Lines_measurements[i-1].Center=P[3*i-2]
		Lines_measurements[i-1].ErrCenter=PCERROR[3*i-2]
		Lines_measurements[i-1].Width=P[3*i-1]
	        Lines_measurements[i-1].ErrWidth=PCERROR[3*i-1]
		Lines_measurements[i-1].Peak=P[3*i]
		Lines_measurements[i-1].ErrPeak=PCERROR[3*i]
		Lines_measurements[i-1].Flux_gaussian=(P[3*i-1]*P[3*i]*SQRT(2*!DPI))
		Lines_measurements[i-1].Err_Flux_gaussian=Error_Flux_GAUSS(P[3*i-2:3*i], PCERROR[3*i-2:3*i])
	   ENDFOR
		Lines_measurements.Background=P[0]
		Lines_measurements.ErrBackground=PCERROR[0]
		Lines_measurements.Chi2=BESTNORM/DOF
  		
  	   FOR i=1,N_gauss do begin
  		;Measure integrated flux under the line
  		Flux_Bkgsub=Flux_line-REPLICATE(Lines_measurements[i-1].Background,n_elements(Flux_line))
		neg_pixel=WHERE(Flux_Bkgsub LT 0,cnt_neg)
		IF cnt_neg GT 1 then  Flux_Bkgsub(neg_pixel)=0. 
  	  	line_region=WHERE(lambda_line GE Lines_measurements[i-1].Center-Lines_measurements[i-1].Width*3.5 $
								and lambda_line LT Lines_measurements[i-1].Center+$
								Lines_measurements[i-1].Width*3.5,N_spread)
  		Lines_measurements[i-1].Flux_integrated=INT_TABULATED(lambda_line(Line_region),Flux_Bkgsub(Line_region))
  		
  		; Compute EW 	
		 	IF KEYWORD_SET(continuum_line) Then begin
		 	    downlimit_region1=max([Lines_measurements[i-1].Center-3.5*Lines_measurements[i-1].Width -default_EW_continuum ,min(lambda_line)])
				uplimit_region1= Lines_measurements[i-1].Center-3.5*Lines_measurements[i-1].Width
				downlimit_region2=Lines_measurements[i-1].Center+3.5*Lines_measurements[i-1].Width
				uplimit_region2=max([Lines_measurements[i-1].Center+ 3.5*Lines_measurements[i-1].Width + default_EW_continuum,max(lambda_line)])    		
				Continnum_region1=WHERE(lambda_line GT downlimit_region1 and $
							lambda_line LT uplimit_region1 ,n_cont1)
				Continnum_region2=WHERE(lambda_line GT downlimit_region2 and $
							lambda_line LT uplimit_region2,n_cont2)
							continuum_region=[Continnum_region1,Continnum_region2]
			    Continuum=mean(continuum_line[continuum_region])
			    Normalized_flux= Flux_line(Line_region)/REPLICATE(Continuum,n_elements(Line_region))
			    Lines_measurements[i-1].EW = INT_TABULATED(lambda_line(Line_region), 1-Normalized_flux)
			 ;   print,Lines_measurements[i-1].Center, Lines_measurements[i-1].EW
			ENDIF
  					
  	   ENDFOR		
  		
	    ;Calculate SNR of detection 
		;calculate level of noise in SN windows
		mean_sigma_line=mean(Lines_measurements.Width)
		downlimit_regionA=max([min(lines_ref)-mean_sigma_line*36,min(lambda_line)])
		uplimit_regionA=min(lines_ref)-mean_sigma_line*3.5
		downlimit_regionB=max(lines_ref)+mean_sigma_line*3.5
		uplimit_regionB=min([max(lines_ref)+mean_sigma_line*36,max(lambda_line)])
		noise_regionA=WHERE(lambda_line GT downlimit_regionA and $
							lambda_line LT uplimit_regionA ,n_noiseA)
						;lambda_line LT uplimit_regionA and error GT 0.,n_noiseA)
		noise_regionB=WHERE(lambda_line GT downlimit_regionB and $
						;lambda_line LT uplimit_regionB and error GT 0.,n_noiseB)
							lambda_line LT uplimit_regionB,n_noiseB)

		
		;Reject lines if the Noise cannot be calculated (region with too many bad pixel)
	    IF n_noiseA LT 2 or n_noiseB LT 2 then begin 		
		Lines_measurements.QC_index = -1
	  	Lines_measurements.Center=-999.
		Lines_measurements.Width=-999.
		Lines_measurements.Peak=-999.
		Lines_measurements.Flux_gaussian=-999.
		Lines_measurements.Flux_integrated=-999.
		Lines_measurements.Flux=-999.
		Lines_measurements.ErrCenter=-999.
	    	Lines_measurements.ErrWidth=-999.
		Lines_measurements.ErrPeak=-999.
		Lines_measurements.Err_Flux_gaussian=-999.

		return,Lines_measurements
	     ENDIF 
	
		sigma_regionA=ROBUST_SIGMA(flux_line(noise_regionA))
		sigma_regionB=ROBUST_SIGMA(flux_line(noise_regionB))
		noise=mean([sigma_regionA,sigma_regionB])
		Lines_measurements.Noise_level=noise
	    	

	  ; Test if the fitted lines are above the detection thresold
	  FOR i=1,N_gauss do begin
		
			SN_detection= Lines_measurements[i-1].Flux_gaussian/(sqrt(N_spread)*noise)
				
				IF SN_detection LT SN_detection_threshold then begin
					Lines_measurements[i-1].QC_index = -3
		  		    Lines_measurements[i-1].Center=-999.
					Lines_measurements[i-1].Width=-999.
					Lines_measurements[i-1].Peak=-999.
					Lines_measurements[i-1].Flux_gaussian=-999.
					Lines_measurements[i-1].Flux=-999.
					Lines_measurements[i-1].Flux_integrated=-999.
				    Lines_measurements[i-1].ErrCenter=-999.
	 			    Lines_measurements[i-1].ErrWidth=-999.
					Lines_measurements[i-1].ErrPeak=-999.
					Lines_measurements[i-1].Err_Flux_gaussian=-999.
					Lines_measurements[i-1].SNR=0.
				    Lines_measurements[i-1].EW=0.
				ENDIF ELSE begin
			 	        Lines_measurements[i-1].SNR=SN_detection
			 		; Test on relative error.
			 		IF Lines_measurements[i-1].Err_Flux_gaussian/Lines_measurements[i-1].Flux_gaussian $
			 								GT MAX_ERROR_FLUX then begin
			  		  Lines_measurements[i-1].QC_index =-2 
			  		  Lines_measurements[i-1].Center=-999.
					  Lines_measurements[i-1].Width=-999.
					  Lines_measurements[i-1].Peak=-999.
					  Lines_measurements[i-1].Flux_gaussian=-999.
					  Lines_measurements[i-1].Flux_integrated=-999.
					  Lines_measurements[i-1].Flux=-999.
				      Lines_measurements[i-1].ErrCenter=-999.
	 			      Lines_measurements[i-1].ErrWidth=-999.
					  Lines_measurements[i-1].ErrPeak=-999.
				      Lines_measurements[i-1].Err_Flux_gaussian=-999.
					  Lines_measurements[i-1].SNR=0.
					  Lines_measurements[i-1].EW=0.

			 		ENDIF ELSE begin
				      Lines_measurements[i-1].QC_index=1 
				      Lines_measurements[i-1].Flux=Lines_measurements[i-1].Flux_integrated
			 		  Lines_measurements[i-1].EW=Lines_measurements[i-1].Flux / Lines_measurements[i-1].Background
			 		ENDELSE
	
			        ENDELSE
	ENDFOR
	FOR i=1,N_gauss do begin
				;Reject line with width over the width limit	
				IF Lines_measurements[i-1].Width GT SIGMA_MAX then begin
					Lines_measurements[i-1].QC_index = -4
		  		    Lines_measurements[i-1].Center=-999.
					Lines_measurements[i-1].Width=-999.
					Lines_measurements[i-1].Peak=-999.
					Lines_measurements[i-1].Flux_gaussian=-999.
					Lines_measurements[i-1].Flux_integrated=-999.
					Lines_measurements[i-1].Flux=-999.
				    Lines_measurements[i-1].ErrCenter=-999.
	 			    Lines_measurements[i-1].ErrWidth=-999.
					Lines_measurements[i-1].ErrPeak=-999.
					Lines_measurements[i-1].Err_Flux_gaussian=-999.
					Lines_measurements[i-1].SNR=0.
					Lines_measurements[i-1].EW=0.

				ENDIF 
				
				;Reject lines affected by strong sky lines
				Line_region=WHERE(lambda_line GE Lines_measurements[i-1].Center-Lines_measurements[i-1].Width*1.5 $
								and lambda_line LT Lines_measurements[i-1].Center+$
								Lines_measurements[i-1].Width*1.5,N_spread)
		                sky_mask=WHERE( mask(line_region) EQ 5, n_sky)
				IF n_sky GE 1 then begin
					Lines_measurements[i-1].QC_index = -5
		  		        Lines_measurements[i-1].Center=-999.
					Lines_measurements[i-1].Width=-999.
					Lines_measurements[i-1].Peak=-999.
					Lines_measurements[i-1].Flux_gaussian=-999.
					Lines_measurements[i-1].Flux_integrated=-999.
					Lines_measurements[i-1].Flux=-999.
				        Lines_measurements[i-1].ErrCenter=-999.
	 			        Lines_measurements[i-1].ErrWidth=-999.
					Lines_measurements[i-1].ErrPeak=-999.
					Lines_measurements[i-1].Err_Flux_gaussian=-999.
					Lines_measurements[i-1].SNR=0.
					Lines_measurements[i-1].EW=0.

				ENDIF 
			    
			   	ENDFOR

	order_lines=sort(lines_ref)
	FOR i=0,N_gauss -2 do begin

				; Test Bended lines
				
				det=WHERE(Lines_measurements.Width GT 0,cnt_det)
				IF cnt_det GT 1 then begin
				mean_width=mean(Lines_measurements.Width(det))
				IF abs(lines_ref[order_lines(i+1)] - lines_ref[order_lines(i)] ) LT 6*mean_width then begin
 				Lines_measurements[order_lines(i)].Flux_integrated=-999.
				Lines_measurements[order_lines(i)].QC_index = 2
 				Lines_measurements[order_lines(i+1)].Flux_integrated=-999.
				Lines_measurements[order_lines(i+1)].QC_index = 2
				Lines_measurements[order_lines(i)].Flux=Lines_measurements[order_lines(i)].Flux_gaussian
				Lines_measurements[order_lines(i+1)].Flux=Lines_measurements[order_lines(i+1)].Flux_gaussian
				
					IF KEYWORD_SET(continuum_line) Then begin
					mean_sigma_line=mean(Lines_measurements.Width)
					downlimit_region1=max([min(lines_ref)-mean_sigma_line*3.5 -default_EW_continuum,min(lambda_line)])
					uplimit_region1=min(lines_ref)-mean_sigma_line*3.5
					downlimit_region2=max(lines_ref)+mean_sigma_line*3.5
					uplimit_region2=min([max(lines_ref)+mean_sigma_line*3.5+default_EW_continuum,max(lambda_line)])
					Continnum_region1=WHERE(lambda_line GT downlimit_region1 and $
								lambda_line LT uplimit_region1 ,n_cont1)
							;lambda_line LT uplimit_regionA and error GT 0.,n_noiseA)
					Continnum_region2=WHERE(lambda_line GT downlimit_region2 and $
							;lambda_line LT uplimit_regionB and error GT 0.,n_noiseB)
								lambda_line LT uplimit_region2,n_cont2)
				  Continuum_region=[Continnum_region1,Continnum_region2]		
				    Continuum=mean(continuum_line[Continuum_region])
				    Normalized_flux= Flux_line(Line_region)/REPLICATE(Continuum,n_elements(Line_region))
				    Lines_measurements[order_lines(i)].EW=Lines_measurements[order_lines(i)].Flux / Continuum
				    Lines_measurements[order_lines(i+1)].EW=Lines_measurements[order_lines(i+1)].Flux / Continuum
				  ;  print,Lines_measurements[order_lines(i)].Center, Lines_measurements[order_lines(i)].EW, '(blended)'
				  ;  print,Lines_measurements[order_lines(i+1)].Center, Lines_measurements[order_lines(i+1)].EW, '(blended)'
					ENDIF
				
				ENDIF
				ENDIF
	ENDFOR

	IF KEYWORD_SET(PLOT_FIT) THEN BEGIN 
	cgPlot,lambda_line, smooth(flux_line ,1)
	cgOplot, lambda_line, error  , color='green'
	cgOplot, lambda_line, gauss_line  , color='blue'
	cgOplot, lambda_line(noise_regionA), noise_regionA*0.  , color='red'
	cgOplot, lambda_line(noise_regionB), noise_regionB*0.  , color='red'
	ENDIF
	
	return,Lines_measurements

	END



	;------------------------------------------------------------------------------------------
	; PROCEDURE  GAMA_FIT_LINES
	;
	; DESCRIPTION : Fit all the emission lines listed in GAMMA_Emission_lines_model.dat in a spectrum
	;				Estimate the 1-sigma flux upper limit for each line, from the level of noise
	; 	
	; INPUT  : Object => Name of the object to fit
	;
	; OUTPUT : Lines_array -> 2D array [fit_group,Center,sigma,flux,errFlux,Noise,SN_detection, FLAG_FIT, PeakLine,Continum_0, Continum_1] x N_lines
	;		   QC_FIT	-> Quality Control flag for emission fit
	;					-> -8 no fitted yet
	;					-> -2 No emission lines fitted for this spectrum
	;					-> -1 No input fits file
	;					->  1 Fit ok
	;-------------------------------------------------------------------------------------------
	
	
	PRO GAMA_FIT_LINES, spectrum_struc,line_library,Measurements_struct,QC_FIT
	
	@GAMA_Main_Config

	;------ Load line libraries  ----------
		lines_name=line_library.Name
		lines_ref=line_library.lines_ref
		lines_group=line_library.lines_group
		
		QC_FIT=-8 ; no fit
		N_lines=n_elements(lines_ref)
		;Make structures for recovering the flux measurement
	 	Measurements_struct=replicate({LineMeasure}, N_lines)
		 	
	;------ Prepare GAMA spectrum ----------
		 		
		 ;Load spectrum 	

		 lambda=spectrum_struc.lambda
		 spectrum=spectrum_struc.spectrum
		 error=spectrum_struc.error
		 mask=spectrum_struc.mask
		 Continuum=spectrum_struc.Continuum
		 object=spectrum_struc.name
		 Spectrum_emission=spectrum-Continuum +replicate(1,n_elements(lambda)) ; add 1 to avoid measurement close to 0
		 Z=SXPAR(spectrum_struc.header,'Z')
		
		
		;Set to 1 data pixels with infinite values in spectrum
		Infinite_spectrum= ~FINITE(Spectrum_emission)
		IF TOTAL(Infinite_spectrum) GE 1 then begin 
			error(WHERE(Infinite_spectrum))=100000;!VALUES.F_NAN
			Spectrum_emission(WHERE(Infinite_spectrum))=1;!VALUES.F_NAN
			mask(WHERE(Infinite_spectrum))=22
		endif
		zero_error= WHERE(error EQ 0., n_zero_error)
		IF  n_zero_error GE 1 then begin 
			error(zero_error)=10000;!VALUES.F_NAN
		endif
		;add error from fitting residual
		 mask_lines=GAMA_LOAD_EMIMASK(lambda)
		 mask_4errors=mask_lines + mask	
		 good_pixel=WHERE(mask_4errors LE 1, n_good)	
		IF  n_good GE 20 then residual=smooth(INTERPOL_V8(Spectrum_emission(good_pixel),lambda(good_pixel),lambda, /NAN),10)	
		IF  n_good LT 20 then residual=lambda * 0.
		 error_total= (sqrt(error^2+residual^2))^2  ;- Normal weighting
		 error_total=error
		
		;Mask regions with strong sky emission
		sky5577=WHERE(lambda GE 5574/(1+Z) and lambda LE 5580/(1+Z), cnt_sky5577) 
		sky6302=WHERE(lambda GE 6299/(1+Z) and lambda LE 6305/(1+Z), cnt_sky6302) 
	        IF cnt_sky5577 GT 1 then mask(sky5577)=5		
  	 	IF cnt_sky6302 GT 1 then mask(sky6302)=5
		

	;------ Fit simultaneously lines in group ----------
			
		FOR i=0, max(lines_group) do begin
		select_lines=WHERE(lines_group EQ i, N_lines_group)
	
		;Test if the line is out off wavelength range
		Lines_out_range=WHERE(lines_ref(select_lines) LT (min(lambda)+20) or $
									lines_ref(select_lines) GT (max(lambda)-20), n_lines_out)
		IF  n_lines_out GE 1 then begin 
		Measurements_struct[Lines_out_range].QC_index=-4
	     	Measurements_struct[Lines_out_range].Name=lines_name[select_lines[Lines_out_range]]
	  
	    	IF n_lines_out LT n_elements(select_lines) then begin
	   		 REMOVE,n_lines_out,select_lines ;Remove from the fit, lines out of range
			ENDIF 
			IF n_lines_out EQ n_elements(select_lines) then begin
			i++
			goto,next_group
			ENDIF 
		ENDIF 
	
			; Trim the spectrum the line window
			region2trim=[max([min(lines_ref(select_lines))-default_Line_window,min(lambda,/NAN)],/NAN), $
					    min([max(lines_ref(select_lines))+default_Line_window,max(lambda,/NAN)],/NAN)]
		 	flux_line=GAMA_TRIM_SPECTRUM(lambda,Spectrum_emission,region2trim)
			lambda_line=GAMA_TRIM_SPECTRUM(lambda,lambda,region2trim)
			error_line=GAMA_TRIM_SPECTRUM(lambda,error_total,region2trim)
			mask_line=GAMA_TRIM_SPECTRUM(lambda,mask,region2trim)
			continuum_line =GAMA_TRIM_SPECTRUM(lambda,Continuum,region2trim)
			
			;Fit doublet
			IF N_lines_group EQ 2 then begin
			line_2_fit=lines_ref(select_lines)
			Measurements_struct[select_lines]=GAMA_multi_gaussian(lambda_line,flux_line,error_line,mask_line,line_2_fit,/FIX_SIGMA,$
											  continuum_line=Continuum)
		        Measurements_struct[select_lines].Name=lines_name[select_lines]
		    
		        ;Fit multi-gaussian or single lines		
			ENDIF ELSE  begin
			line_2_fit=lines_ref(select_lines)
			Measurements_struct[select_lines]=GAMA_multi_gaussian(lambda_line,flux_line,error_line,mask_line,line_2_fit, $
											  continuum_line=Continuum)
			Measurements_struct[select_lines].Name=lines_name[select_lines]
			ENDELSE
			Measurements_struct[select_lines].Background=Measurements_struct[select_lines].Background -1 ; subtract -1 added to the input spectrum
			
next_group:			
		
		ENDFOR
	
	   QC_FIT=1 ;Sucessfully fit all the groups 
	
	   ;test if all the fits had failed
	   failure_array=Measurements_struct[*].QC_index
	   tmp=WHERE(failure_array LT 0, n_lines_failure)
	   IF n_lines_failure EQ N_lines then begin
	   QC_FIT=-2
	   return
	   ENDIF
	
	   ; Calculate upper flux limit for undetected lines
	   Detected_lines=WHERE(Measurements_struct[*].SNR GE SN_detection_threshold)
	   Undetected_lines=WHERE(Measurements_struct[*].SNR LT SN_detection_threshold, N_Undetected_lines) 
	   IF N_Undetected_lines EQ 1 then begin
		Upper_limit_flux=sqrt(2*pi)*Measurements_struct[N_Undetected_lines].Noise_level*SIGMA_DEFAULT
		Measurements_struct[Undetected_lines].Flux_limit=Upper_limit_flux
		;Measurements_struct[Undetected_lines].SNR=0.

	    ENDIF
	 
	   IF N_Undetected_lines GE 2 then begin
		Mean_sigma_all=mean(Measurements_struct[Detected_lines].Width)
		Upper_limit_flux=sqrt(2*pi)*Measurements_struct[Undetected_lines].Noise_level*replicate(Mean_sigma_all,N_Undetected_lines)
		Measurements_struct[Undetected_lines].Flux_limit=Upper_limit_flux[Undetected_lines]
		;Measurements_struct[Undetected_lines].SNR=0.
	   ENDIF

	return
	END


	;------------------------------------------------------------------------------------------
	; ROUTINE  GAMA_FitEmi_all
	;
	; DESCRIPTION : Main routine to measure the emission lines of all GAMA spectra
	; 	
	; INPUT  : none
	;
	; OUTPUT : Save file, with 
	;          Array_line_all -> 2D array N_lines x [fit_group,Center,sigma,flux,errFlux,SN_detection, FLAG_FIT] x N_target
	;		   QC_FIT	-> Quality Control flag for emission fit
	;					-> -9 no fitted yet
	;					-> -2 No emission lines fitted for this spectrum
	;					-> -1 No input fits file
	;					->  1 Fit ok
	;-------------------------------------------------------------------------------------------
	
	
	
 PRO GAMA_FitEmi_all
	
 @GAMA_Main_Config
	
 print,"Saved in:"+PATH_CATALOGUE+Saved_Lines
 readcol,QC1_CATALOGUE_OUTPUT+'.dat',format='(a,f,f,a,i,a)', SPECID,CATAID,z,QC_Flag,QF,Q1_Flag
 good_to_fit=WHERE(STRCMP(Q1_Flag,'G'),N_target)
 targets=SPECID(good_to_fit)
 QC2_FIT=intarr(N_target)
 QC2_FIT(*)=-9
	
 ; Load library of lines to fit 
 IF (FILE_EXIST('GAMMA_Emission_lines_model.dat')) then begin
		readcol,"GAMMA_Emission_lines_model.dat",format='(a,f,i)',lines_name, lines_ref,lines_group,/silent
		Line_library={Name:lines_name, lines_ref:lines_ref,lines_group:lines_group}	
 ENDIF ELSE begin 
 print, "Catalogue of index lines not found " 
 stop
 ENDELSE 

	
 ; Create a structure to save all the measurements
 N_lines=n_elements(lines_ref)
 Line_all=replicate({LineMeasure}, N_lines,1)

	
 timer_start= SYSTIME(1)  
	
 FOR i=0, N_target - 1 do begin
		
	object=targets(i)

	file=preffixe_final+object+".fits"
	 
	Spec_struc=GAMA_READ_FITS(PATH_FIT+file,/SILENT)
	check_data=size(Spec_struc) ;Check the data size 
		IF (check_data(0) EQ 0) then begin QC2_FIT(i)=-1 ; no file
			print,'Failed to open file:',PATH_FIT+file
		ENDIF ELSE begin 
			GAMA_FIT_LINES, Spec_struc,Line_library,Measurements_struct,QC_FIT_target
			Line_all=[[Line_all],[Measurements_struct]]
			QC2_FIT[i]=QC_FIT_target	
		ENDELSE
		
	;Autosave
	IF ((i MOD AUTO_SAVE_Lines) EQ 0) then begin
	save, targets,Line_all, QC2_FIT, filename=PATH_CATALOGUE+Saved_Lines
	ENDIF
		
	;print progression
	IF ((i MOD 100) EQ 0) then begin
	tmp=WHERE(QC2_FIT EQ 1,N_Good_fit)
	tmp=WHERE(QC2_FIT GT -9 and QC2_FIT LE -1,N_Bad_fit)
	print,"["+STRTRIM(i,2)+"/"+STRTRIM(N_target,2)+"] :"+STRTRIM(N_Good_fit,2)+" succesfull | "+STRTRIM(N_Bad_fit,2)+" Failure"
	ENDIF
		
 ENDFOR

 save, targets,Line_all, QC2_FIT, filename=PATH_CATALOGUE+Saved_Lines
	
 print, "   [FIT FINISH] "+STRTRIM(ceil(SYSTIME(1)-timer_start),2)+"sec"
 print, STRTRIM(ceil(SYSTIME(1)-timer_start)/N_target,2)+"sec / target"
	
	
 END
