	;-----------------------------------------------------------------------
	; PROCEDURE Plot_all_lines
	; DESCRIPTION :  Plot all the lines and the fit in a multiplot 
	;				 Saved in a .png file
	; INPUT		  : Lines_array -> 2D array [fit_group,Center,sigma,flux,errFlux,Noise $
	;                 		,SN_detection, FLAG_FIT, PeakLine,Continum_0, Continum_1] x N_lines
	;				Object -> String containing the name of the object
	;				Lambda -> wavelength array
	;				Flux -> Flux array
	;				Weight -> Weight array
	; 				
	; OUTPUT : NONE
	;---------------------------------------------------------------------
	
	PRO Plot_all_lines,Object,Lines_array, lambda,Spectrum,error, mask
	
	@GAMA_Main_Config
	
	readcol,"GAMMA_Emission_lines_model.dat",format='(a,f,i)',lines_name, lines_ref,lines_group,/silent
		;print, "Group  Line  |   Lambda   Sigma   Flux   SN_detec"
	
	N_group=max(lines_group)
	
	;PS_Start,PATH_LINES+'Fig/'+Object+'_lines.ps',/quiet
	;cgDisplay, 700,600
	positions = cgLayout([3,3], OXMargin=[3,3], OYMargin=[2,4], XGap=4, YGap=4)
	;positions = cglayout([3,3], XGap=4, YGap=4)
	
	cgPS_Open, PATH_LINES+'Lines_'+Object+'.eps'
	FOR i=0, N_group do begin
	
		Select_Lines=WHERE(lines_group EQ i, N_lines_group)
		Continuum=Lines_array[Select_Lines].Background 
		Center_ref=lines_ref[Select_Lines]
		Center_gauss=Lines_array[Select_Lines].Center 
		Sigma_gauss=Lines_array[Select_Lines].Width
		Peak_gauss=Lines_array[Select_Lines].Peak  
		SN_detection=Lines_array[Select_Lines].SNR  
		default_window=default_Line_window*0.8
		region2trim=[max([min(Center_ref)-default_window,min(lambda,/NAN)])$
	 				,min([max(Center_ref)+default_window,max(lambda,/NAN)])]


	 flux_line=GAMA_TRIM_SPECTRUM(lambda,Spectrum,region2trim)
	 lambda_line=GAMA_TRIM_SPECTRUM(lambda,lambda,region2trim)
	 error_line=GAMA_TRIM_SPECTRUM(lambda,error,region2trim)
	 mask_line=GAMA_TRIM_SPECTRUM(lambda,mask,region2trim)
	

		 fit_line=lambda_line*0+Continuum[0]

		FOR j=0, N_lines_group-1 do begin
		print,Lines_array[Select_Lines[j]].Name ,Center_gauss[j], Sigma_gauss[j], SN_detection[j],Lines_array[Select_Lines[j]].QC_index,Continuum[j]
			IF Center_gauss[j] NE -999.000 then begin
			fit_line= fit_line+ GAUSS1(lambda_line, [Center_gauss[j], Sigma_gauss[j], Peak_gauss[j]],/PEAK)
			
			endif
		ENDFOR
	
	       ;Set to 1 data pixels with infinite values in spectrum
		Infinite_spectrum= ~FINITE(flux_line)
		IF TOTAL(Infinite_spectrum) GE 1 then begin 
			error_line(WHERE(Infinite_spectrum))=100000;!VALUES.F_NAN
			flux_line(WHERE(Infinite_spectrum))=1;!VALUES.F_NAN
			mask(WHERE(Infinite_spectrum))=22
		endif
		zero_error= WHERE(error_line EQ 0., n_zero_error)
		IF  n_zero_error GE 1 then begin 
			error_line(zero_error)=10000;!VALUES.F_NAN
		endif
		;add error from fitting residual
		 mask_lines=GAMA_LOAD_EMIMASK(lambda_line)
		 mask_4errors=mask_lines + mask_line	
		 good_pixel=WHERE(mask_4errors LE 1, n_good)	
		IF n_good GT 10 then begin
		 residual=smooth(INTERPOL_V8(flux_line(good_pixel),lambda_line(good_pixel),lambda_line, /NAN),10)	
		 error_total= (sqrt(error_line^2+residual^2))^2  ;- Normal weighting
		endif else begin
		error_total=error_line
		endelse
	;

	;
	pixel_sky=WHERE(mask_line lT 5., cnt_sky)
	IF cnt_sky GT 1 then mask_line(pixel_sky)=0.
	yrange_=[mean(flux_line,/NAN)-2*sqrt(VARIANCE(flux_line,/NAN)),max(flux_line,/NAN)]
	oplotObj1 = Obj_New('cgOverPlot',lambda_line, error_line  , color='green') 
	oplotObj3 = Obj_New('cgOverPlot',lambda_line, mask_line  , color='blue') 
	oplotObj2 = Obj_New('cgOverPlot',lambda_line, fit_line  , color='red',thick=3) 
	cgPlot,lambda_line, flux_line ,OPLOTS=[oplotObj1,oplotObj2,oplotObj3],xrange=region2trim,xs=1, Position=positions[*,i], NoErase=(i NE 0),title=string(lines_name(Select_Lines[0])),charsize=0.9,yrange=yrange_ , XTICKINTERVAL=50, xthi=1, ythi=1, thick=1.5, charthic=1
	endfor
	cgText, 0.5, 0.98, Alignment=0.5, /Normal, Object, color='royal blue',charsiz=1.4, charthic=1
cgPS_Close, /PNG,  Width=900
	return

	END
	
	
	
	PRO GAMA_Plot_All_Objects
	
	@GAMA_Main_Config

	restore,PATH_CATALOGUE+saved_lines
	readcol,QC1_CATALOGUE_OUTPUT+'.dat',format='(a,f,f,a,i,a)', SPECID,CATAID,z,QC_Flag,QF,Q1_Flag

	good_to_fit=WHERE(STRCMP(Q1_Flag,'G'),N_target)
	targets=SPECID(good_to_fit)

	for i=0, n_elements(targets) -1 do begin
	file=preffixe_final+targets(i)+".fits"
print,PATH_FIT+file
	Spec_struc=GAMA_READ_FITS(PATH_FIT+file,/SILENT)
	lambda=Spec_struc.lambda
	Spectrum=Spec_struc.Spectrum - Spec_struc.Continuum
	error=Spec_struc.error
	mask=Spec_struc.mask
	Lines_array=LINE_ALL[*,i+1]
	
	;Mask regions with strong sky emission
	sky5577=WHERE(lambda GE 5570/(1+Z(i)) and lambda LE 5582/(1+Z(i)), cnt_sky5577) 
	sky6302=WHERE(lambda GE 6298/(1+Z(i)) and lambda LE 6306/(1+Z(i)), cnt_sky6302) 
	 IF cnt_sky5577 GT 1 then mask(sky5577)=5		
  	 IF cnt_sky6302 GT 1 then mask(sky6302)=5

	Plot_all_lines,targets(i),Lines_array, lambda,Spectrum,error, mask
;PAUSE
	endfor
	
	END
