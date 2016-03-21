pro IndexMeasure__define
 
;  This routine defines the structure for the index measurement grid
 
  tmp = {IndexMeasure, $
         Name: ' ',$         		 ; Index Name
         Index: -999.d,   $     	 ; Index Value
         Index_err: -999.d, $   	 ; Index err
         Continum_bleu: -999.d,   $  ; Continum_bleu
         Continum_bleu_error: -999.d , $  ; Continum bleu error
         Continum_red:   -999.d,   $      ; Continum_red
         Continum_red_error: -999.d,   $  ; Continum red error
         Type: 'x',    $   				  ; Type of index (index, with emission, break)
         Model_Flag: -9, $        		  ; Flag model spectrum (0 = No model, 1 = Model)
         QC_index:-99 $ 				  ; Quality Flag
         }
 
end


PRO GAMA_DOWNGRADE_RESOLUTION, Spectrum
@GAMA_Main_Config

            ;-----------------------------------------------------------
            ; Find the size for the gaussian -- we want a gaussian that
            ; extends to 4.0 sigma (99.99%), but it has to have an odd
            ; number of pixels in order to be properly centered.  
            ; Force the gaussian filter size to be 8*sigma rounded to 
            ; the nearest odd number.
            ; FROM EZ_LICK soft
            ;-----------------------------------------------------------

				vel_disp=0.
				
                dlambda = (max(Spectrum.Lambda)-min(Spectrum.Lambda))/((size(Spectrum.Lambda))[1]-1)
                lambda = mean(Spectrum.Lambda)

                sig_ind = LICK_RESOLUTION / sqrt(8.0*alog(2.0)) / dlambda
                sig_gal = (FLOAT(vel_disp) / c) * (lambda / dlambda)
       
                sig_tot = sqrt(NOMINAL_RESOLUTION^2 + sig_gal^2)
                
                if (sig_ind GT sig_tot) then begin

                  sigma = sqrt(sig_ind^2 - sig_tot^2)
                  filtersize=(FIX(sigma*8)+(FIX(sigma*8+1) mod 2))
                  gaussfilter = PSF_GAUSSIAN(NPIXEL=filtersize, $
                                         ST_DEV=sigma, NDIMEN=1,/NORMALIZE)
                  smoothspec = CONVOL(Spectrum.spectrum, gaussfilter, /CENTER, $
                                        /EDGE_TRUNCATE)
	
                 smooth_err = CONVOL((Spectrum.Error)^2, gaussfilter, /CENTER, $ ;Smoothing reduces error 
                                        /EDGE_TRUNCATE) 
                                        ; add gaussian-weighted errors from surrounding pixels in quadrature
                 smooth_err = SQRT(smooth_err) 
                                     
				Spectrum.spectrum=smoothspec
				Spectrum.Error=smooth_err
					
				endif

END



FUNCTION GAMA_Measure_index, Spectrum,index_BP,MODEL=MODEL,IND_BREAK=IND_BREAK,EMISSION=EMISSION, $  
												magnitude=magnitude,PLOT_LINE=PLOT_LINE

@GAMA_Main_Config

  ;Initialize the measure structure
    
      mesure_index={IndexMeasure}
	  mesure_index.Name=index_BP.name
      mesure_index.Index=-999
      mesure_index.Index_err=-999
      mesure_index.Continum_bleu=-999
      mesure_index.Continum_red=-999
      mesure_index.Continum_bleu_error=-999
      mesure_index.Continum_red_error=-999
      mesure_index.Type='x'
	  mesure_index.Model_Flag='9'
  	  mesure_index.QC_index=-99; QC flag

  
  IF KEYWORD_SET(MODEL) THEN BEGIN ;Fitting stellar model
  	 IF ~(TAG_EXIST(Spectrum, 'Stellar_continum_model', /QUIET  )) then begin 
   		print, "/MODEL keyword set but the input fit file does not have a stellar continuum model array"
   		stop
   		mesure_index.QC_index=-9
   		return, mesure_index
   	  ENDIF
   	
  Flux=Spectrum.Stellar_continum_model
  ErrFlux=1/sqrt(Spectrum.Stellar_continum_model)
  mask=fltarr(n_elements(Flux))
  mask(*)=0 
  good=WHERE(Spectrum.mask EQ 0)
  Model_Flag=1

  ENDIF ELSE BEGIN ;Fitting a raw data
  Flux=Spectrum.Spectrum
  ErrFlux=Spectrum.error
  good=WHERE(Spectrum.mask EQ 0)
  mask=Spectrum.mask 
  Model_Flag=0
  ENDELSE
  
; Test bad pixels
  	;Verify that the index is not out of the spectrum range 
	if (index_BP.Blue[0] LT min(Spectrum.lambda)) and (index_BP.Red[1] GT max(Spectrum.lambda)) then begin
	IF PRINT_RESULTS THEN print,"Index line out of the spectral range"
	mesure_index.QC_index=-8
	return, mesure_index
	endif
  
    ;Test the amount of bad pixel in the line region
	Lick_region=WHERE(Spectrum.lambda GE index_BP.Blue[0] $
	 					and Spectrum.lambda LE index_BP.Red[1],N_elements_lines)
	tmp=WHERE(mask[Lick_region] EQ 0.,N_Good_line )
	IF (N_Good_line+1.)/N_elements_lines LT (1-BP_MAX_FRACTION) then begin  
	mesure_index.QC_index=-5
	IF PRINT_RESULTS THEN print,"Too many wrong pixel in the region"
	return, mesure_index
	ENDIF
     

 ; Measure continuum level
 Type='i'
      Continum_bleu=im_integral(Spectrum.lambda[good],Flux[good],$
      						index_BP.Blue[0],index_BP.Blue[1])/index_BP.Blue[2]
      Continum_bleu_error=sqrt(im_integral(Spectrum.lambda[good],ErrFlux[good]^2,$
      						index_BP.Blue[0],index_BP.Blue[1]))/index_BP.Blue[2]	
    
      Continum_red=im_integral(Spectrum.lambda[good],Flux[good],index_BP.Red[0]$
      				,index_BP.Red[1])/index_BP.Red[2]
      Continum_red_error=sqrt(im_integral(Spectrum.lambda[good],ErrFlux[good]^2, $ 
      					index_BP.Red[0],index_BP.Red[1]))/index_BP.Red[2]   
   
    BPc_center = [index_BP.Blue[0] + index_BP.Blue[2]/2.0,index_BP.Red[0] + index_BP.Red[2]/2.0]

    ; if both continuum are defined, fit a polynomial between the two
    ; pseudo-continuum points; to get the errors in the coefficients right
    ; subtract off the starting wavelength
	
	;Calculate index for Breaks
 	if keyword_set(IND_BREAK)  then begin
 		break_index=(index_BP.Blue[2]/index_BP.Red[2])*(Continum_red/Continum_bleu)
 		err_break_index=0.
 		IF PRINT_RESULTS THEN print,index_BP.name,' | ',TRIM(break_index,'(f6.2)',2)
 		index=break_index
 		index_err=err_break_index
 		Type='b'
 		GOTO, Plot_and_RETURN
 	endif

	;Test the amount of bad pixel in the line region
	Line_region=WHERE(Spectrum.lambda GE index_BP.Line[0] $
	 					and Spectrum.lambda LE index_BP.Line[1],N_elements_lines)
	tmp=WHERE(mask[Line_region] EQ 0,N_Good_line )
	IF (N_Good_line+1.)/N_elements_lines LT (1-BP_MAX_FRACTION) then begin  
	mesure_index.QC_index=-5
	IF PRINT_RESULTS THEN print,"Too many wrong pixel in the line region"
	return, mesure_index
	ENDIF
      
	;Compute continuum using a least square fiting
    if (Continum_bleu_error gt 0.0) and (Continum_red_error gt 0.0) then begin
       Continuum_coef = poly_fit(BPc_center,[Continum_bleu,Continum_red],$
        1,measure_errors=[Continum_bleu_error,Continum_red_error],sigma=ccoeff_err,$
         /double,chisq=chisq,covar=covar,yerror=yerror,yband=yband,$
         yfit=yfit,status=status)
        Continuum_err_coef=ccoeff_err
    endif else begin  
    	 mesure_index.QC_index=-7 ; Fit of the continuum failed
     	 return, mesure_index
    endelse
	
     ;Subtract the emission line 
	 IF KEYWORD_SET(EMISSION) and ~(KEYWORD_SET(MODEL)) THEN BEGIN 
    	 Spectrum_emission=Flux-Spectrum.Stellar_continum_model
   	     emi_line=GAUSSFIT(Spectrum.lambda(Line_region),Spectrum_emission(Line_region),Gauss_param,NTERMS=3)
    	 Flux(Line_region)=Spectrum.spectrum(Line_region)-emi_line
     	 Type='e'
     ENDIF
     		
	; compute the line EW; evaluate the continuum over the absorption line
	; and continuum points; see Trager et al. (1998)

    if (Continum_bleu_error gt 0.0) and $
      (Continum_red_error gt 0.0) and (index_BP.Line[0] ge min(Spectrum.lambda)) and $
      (index_BP.Line[1] le max(Spectrum.lambda)) then begin
      
       ;Compute Pseudo continum and measure flux under the line
       clineflux = poly(Spectrum.lambda,Continuum_coef)
       clineflux_err = sqrt(poly((Spectrum.lambda)^2,Continuum_err_coef^2)) ; continuum error
       cinline=im_integral(Spectrum.lambda[good],(clineflux)[good],$
       			index_BP.Line[0],index_BP.Line[1])/(index_BP.Line[1]-index_BP.Line[0])
       cinline_err = sqrt(im_integral(Spectrum.lambda[good],clineflux_err[good]^2, $ 
               index_BP.Line[0],index_BP.Line[1]))/(index_BP.Line[1]-index_BP.Line[0])

       ;Flux in the line
       lineflux = im_integral(Spectrum.lambda[good],(clineflux-Flux)[good],$
       			index_BP.Line[0],index_BP.Line[1])/(index_BP.Line[1]-index_BP.Line[0])
       lineflux_err = sqrt(im_integral(Spectrum.lambda[good],ErrFlux[good]^2,$
       				index_BP.Line[0],index_BP.Line[1]))/(index_BP.Line[1]-index_BP.Line[0])
         
       ; EW of the line 
       linefluxnorm = 1.0 - Flux/clineflux
       linefluxnorm_err = ErrFlux/clineflux
          
       index = im_integral(Spectrum.lambda[good],linefluxnorm[good],index_BP.Line[0],index_BP.Line[1])
       index_err = sqrt(im_integral(Spectrum.lambda[good],linefluxnorm_err[good]^2,$
       			   index_BP.Line[0],index_BP.Line[1]))
   	  
   	   IF PRINT_RESULTS THEN  print,index_BP.name,' | ',TRIM(index,'(f6.2)',2),'+/-',TRIM(index_err,'(f6.2)', $
   	   		2),' | ', TRIM(lineflux,'(f6.2)',2),'+/-',TRIM(lineflux_err,'(f6.2)',2),' | ',$
   	   		TRIM(Continum_bleu,'(f6.2)',2),'+/-',TRIM(Continum_bleu_error,'(f6.2)',2),' | ', $
   	   		TRIM(Continum_red,'(f6.2)',2),'+/-',TRIM(Continum_red_error,'(f6.2)',2)

    endif else begin  
    mesure_index.QC_index=-6
    return, mesure_index
    endelse

    
Plot_and_RETURN:
    
  IF KEYWORD_SET(PLOT_LINE) THEN BEGIN 
     cgPS_Open,path_figure+Spectrum.name+'_'+index_BP.name+'.ps'

    IF ~(KEYWORD_SET(MODEL)) then Flux=smooth(Flux,SMOOTH_4PLOT)
    define_xrange=[index_BP.Blue[0]-100,index_BP.Red[1]+100]
    define_yrange=[min(Flux(WHERE(Spectrum.lambda GT define_xrange[0] and $
    				Spectrum.lambda LT define_xrange[1])),/NAN)*0.8$
    			   ,max(Flux(WHERE(Spectrum.lambda GT define_xrange[0] $
    			   and Spectrum.lambda LT define_xrange[1])),/NAN)*1.2]
     
     position = [0.125, 0.125, 0.9, 0.9]
     transparent = 70 ; Percent of transparency.
     

	 thisDevice = !D.Name				 
     Set_Plot, 'Z'
     Device, Set_Resolution=[1240,1010], Decomposed=1, Set_Pixel_Depth=24
     cgErase
     cgplot,Spectrum.lambda,Flux,xrange=define_xrange,yrange=define_yrange,XS=4,YS=4, Position=position
     Blue_region_plot=[[index_BP.Blue[0],index_BP.Blue[1],index_BP.Blue[1],index_BP.Blue[0],index_BP.Blue[0]],$
					 [!y.crange[0],!y.crange[0],!y.crange[1],!y.crange[1],!y.crange[0]]]
	 Red_region_plot=[[index_BP.Red[0],index_BP.Red[1],index_BP.Red[1],index_BP.Red[0],index_BP.Red[0]],$
					 [!y.crange[0],!y.crange[0],!y.crange[1],!y.crange[1],!y.crange[0]]]
     cgoplot,BPc_center,[Continum_bleu,Continum_red], SymSize=2, PSym=cgSYMCAT(45),color='violet', $
     	    Position=position	
     	    
     	IF ~(KEYWORD_SET(IND_BREAK)) THEN BEGIN
     		cgoplot,Spectrum.lambda,clineflux,color='violet'   , Position=position
	    ENDIF
	    
     background_image = cgSnapshot(Position=position)
     cgColorFill, [Spectrum.lambda,REVERSE(Spectrum.lambda)],[(!y.crange[1])*mask,0*mask], color='grey'
   	 cgColorFill, Blue_region_plot[*,0], Blue_region_plot[*,1], Color='blu6'
   	 cgColorFill, Red_region_plot[*,0], Red_region_plot[*,1], Color='Rose'
     foreground_image = cgSnapshot(Position=position)
  	 Set_Plot, thisDevice 
     cgDisplay, 1240,1010 
     alpha = 1 - (transparent/100.)
     cgBlendimage, foreground_image, background_image, ALPHA=alpha, POSITION=position 	
   	 cgplot,Spectrum.lambda,Flux,xrange=define_xrange,yrange=define_yrange,XS=1,YS=1,xtitle='Wavelenght [A]',ytitle='Normalized flux',title=index_BP.name,/NoData, Position=position , /NoErase

	cgPS_Close
		
  ENDIF 

       if keyword_set(magnitude) and (index gt 0.0) then begin
          index_err = (2.5/alog(10.0))*index_err/index
          index = -2.5*alog10(index)
       endif
      mesure_index.Index=index
      mesure_index.Index_err=index_err
      mesure_index.Continum_bleu=Continum_bleu
      mesure_index.Continum_red=Continum_red
      mesure_index.Continum_bleu_error=Continum_bleu_error
      mesure_index.Continum_red_error=Continum_red_error
      mesure_index.Type=Type
	  mesure_index.Model_Flag=Model_Flag
	  mesure_index.QC_index=1
      return, mesure_index
END



PRO GAMA_MEASURE_LIST_INDEX,Spectrum, Measurements_struct;,SILENT=SILENT

@GAMA_Main_Config

  ;Load the index config file
	IF (FILE_EXIST('GAMA_index_model.dat')) then begin
	readcol, 'GAMA_index_model.dat',format='(a,a,f,f,f,f,f,f)',$
 				Index,code,B_BP_low,B_BP_high,C_BP_low,C_BP_high,R_BP_low,R_BP_high,/SILENT
	ENDIF ELSE begin 
	print, "Catalogue of index lines not found " 
	stop
	ENDELSE 

	;Make structures for recovering the index measurement
	N_index=n_elements(Index)
 	Measurements_struct=replicate({IndexMeasure}, N_index)

	; Reduce to Lick resolution if requiere
	IF LICK_INDEX THEN GAMA_DOWNGRADE_RESOLUTION, Spectrum
	
    IF PRINT_RESULTS THEN print,"INDEX |  EW [A]    | Flux_line  | Continum_B | Continum_R "
   	   		 	   		
	FOR i=0, N_index-1 do begin 

	index_BP={Name:Index[i], Blue:[B_BP_low[i],B_BP_high[i],(B_BP_high[i]-B_BP_low[i])],Line:[C_BP_low[i],C_BP_high[i],(C_BP_high[i]-C_BP_low[i])],Red:[R_BP_low[i],R_BP_high[i],(R_BP_high[i]-R_BP_low[i])]}

	case code[i] of 
	'i' : begin 			
			index_value=GAMA_Measure_index(Spectrum,index_BP,PLOT_LINE=PLOT_SAVE,MODEL=FLAG_FIT_MODEL) ;
			end
    'b' : begin 			
			index_value=GAMA_Measure_index(Spectrum,index_BP,/IND_BREAK,PLOT_LINE=PLOT_SAVE, MODEL=FLAG_FIT_MODEL) ;,
			end
	'e' : begin 			
			index_value=GAMA_Measure_index(Spectrum,index_BP,PLOT_LINE=PLOT_SAVE,MODEL=FLAG_FIT_MODEL);EMISSION=FLAG_FIT_EMI,
			end
	endcase
	Measurements_struct[i]=	index_value
	ENDFOR	
			
END



PRO GAMA_FIT_ALL

@GAMA_Main_Config

	IF (FILE_EXIST('GAMA_index_model.dat')) then begin
	readcol, 'GAMA_index_model.dat',format='(a,a,f,f,f,f,f,f)',$
 				Index,code,B_BP_low,B_BP_high,C_BP_low,C_BP_high,R_BP_low,R_BP_high,/SILENT
	ENDIF ELSE begin 
	print, "Catalogue of index lines not found " 
	stop
	ENDELSE 

	;Catalogue='list'
	;readcol,path+Catalogue,format='(a)',targets
	
	@GAMA_Main_Config
	
	readcol,PATH_CATALOGUE+CATALOGUE_INPUT_FIT, format="a,l,f,a",spectrumid, cataid, z, QC_Flag
	n_obj=n_elements(spectrumid)
	QF=intarr(n_obj)
	QF(*)=-9 ; not fitted yet

	TO_FIT=WHERE( STRCMP(QC_Flag, 'G') EQ 1, N_obj) 


	;N_obj=n_elements(spectrumid)
	N_lines=n_elements(Index)
	
	print, "--------------------------------------------- "
	print, "-        GAMA Lick Index routine       -"
	print, "  "
	print, " "
	print, " Data DIR : "+PATH_OBS
	print, "  Line catalogue :"+"GAMA_index_model.dat"
	print, "  Catalogue loaded :"+CATALOGUE_INPUT_FIT
	print, "  Number of index to fit: "+TRIM(N_lines,2)
	print, "  Number of spectrum: "+TRIM(N_obj,2)
	print, "--------------------------------------------- "
	print, " "
	print, " "


	
	Results_array=strarr(N_lines*3+1,N_obj)
	Header=strarr(N_lines*3+1)
	Header[0]='"Name"'

FOR jj=0, N_obj-1 do begin

	file=spectrumid(TO_FIT(jj))+".fits"
	Results_array[0,jj]=spectrumid(jj)
	Spec_struc_input=GAMA_READ_FITS2(PATH_OBS+file)
	check_data=size(Spec_struc_input)
	
	IF (check_data(0) EQ 0) then begin 
	print, "        !! Failed to open file "+file
	QF(TO_FIT(jj))=0	
	ENDIF ELSE begin 

	 Spec_struc=GAMA_Prepare_Lick( Spec_struc_input, /TRIM, /MASK_SLICER)
	 
     GAMA_MEASURE_LIST_INDEX,Spec_struc,Measurements_struct
 

    	for k=1, N_lines do begin
   	 	Results_array[3*k-2,jj]=TRIM(Measurements_struct[k-1].Index,2)
    	Results_array[3*k-1,jj]=TRIM(Measurements_struct[k-1].Index_err,2)
    	Results_array[3*k,jj]=TRIM(Measurements_struct[k-1].QC_index,2)
    	IF jj EQ 0 then begin
    	Header[3*k-2]='"Index_'+Measurements_struct[k-1].Name+'"'
    	Header[3*k-1]='"ErrIndex_'+Measurements_struct[k-1].Name+'"'
    	Header[3*k]='"QC_'+Measurements_struct[k-1].Name+'"'
    	ENDIF
    	endfor

	ENDELSE	
				

	IF ((jj MOD SAVE_2_CATALOGUE) EQ 0) then begin
 	save, Header,Results_array, filename=PATH_CATALOGUE+saved_file_index
	ENDIF
	

    ;print progression
	IF ((jj MOD 100) EQ 0) then begin
 	print,"["+STRTRIM(jj,2)+"/"+STRTRIM(N_obj,2)+"] "
	ENDIF

	
ENDFOR

 IF SAVE_2_CATALOGUE THEN BEGIN 
 restore, PATH_CATALOGUE+saved_file_index
 Write_CSV_Data, Results_array, Header,filename= PATH_CATALOGUE+Index_catalogue 
 print,"Index measurements save in "+PATH_CATALOGUE+Index_catalogue
 ENDIF

END