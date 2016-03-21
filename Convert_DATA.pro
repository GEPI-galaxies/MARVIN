; Convert SDSS spectra file into GAMA spectra format
;
;
FUNCTION GAMA_Uniform_resolution, Spec_struc


Blue_arm=[3750.,5650.]
slicer=[5650.,5900.]


	smooth_spectrum=Spec_struc.spectrum
	
	fsmooth_blue = gaussfold(Spec_struc.lambda,Spec_struc.spectrum,1.24, LAMMIN=min(Spec_struc.lambda,/NAN),LAMMAX=Blue_arm[1])
	fsmooth_slicer= gaussfold(Spec_struc.lambda,Spec_struc.spectrum,1.12, LAMMIN=slicer[0],LAMMAX=slicer[1])
	 
	smooth_spectrum(WHERE(Spec_struc.lambda LT Blue_arm[1]))= fsmooth_blue(WHERE(Spec_struc.lambda LT Blue_arm[1]))
	smooth_spectrum(WHERE(Spec_struc.lambda GE slicer[0] and Spec_struc.lambda LT slicer[1]))= fsmooth_slicer(WHERE(Spec_struc.lambda GE slicer[0] and Spec_struc.lambda LT slicer[1]))
	
	;oplotObj = Obj_New('cgOverPlot', Spec_struc.lambda, smooth_spectrum, color='red') 
	;cgPlot, spec_struc.lambda,Spec_struc.spectrum, OPLOTS=oplotObj 
	
	;header=Spec_struc.Header
	;sxaddpar,header,'DLa','1.6'

	Spec_struc.spectrum=smooth_spectrum
	;Spec_struc.Header=header
	
	return, Spec_struc

END

FUNCTION Convert_SDSS_2_GAMA, file , SPECID, FILE_OUT,SILENT=SILENT

IF (FILE_EXIST(file)) then begin
	Data=readfits(file,header_origin, EXTEN_NO=0, /SILENT)	
	 	;Check the data size 
	 	check_data=size(DATA)
			IF (check_data(0) NE 2) then begin 
			 	IF ~(KEYWORD_SET(SILENT)) then begin 
			 	 print,"Wrong fits format"
			     return,0
		    	ENDIF
		    ENDIF	
ENDIF ELSE begin 
		IF ~(KEYWORD_SET(SILENT)) then begin 
		   print,"File not found : "+file
		   return,0
        ENDIF
ENDELSE
	
  OBJECT_NAME= SPECID
  spectrum_ln=Data[*,0]
  error_ln=Data[*,2]
  Mask_ln=Data[*,3]
  
  CRVAL=SXPAR(header_origin,'CRVAL1')
  NAXIS1=SXPAR(header_origin,'NAXIS1')
  CDELT=SXPAR(header_origin,'CD1_1')
  CRPIX=SXPAR(header_origin,'CRPIX1')
  lambda=10^(CRVAL+CDELT*indgen(NAXIS1))
  OBJECT_NAME= SPECID ;SXPAR(header_origin,'SPECID')	

  ;Resample Wavelength in linear
   NAXIS1_new=NAXIS1*2
   cdelt_new=(max(lambda)-min(lambda))/NAXIS1_new
   lambda_linear=indgen(NAXIS1_new)*cdelt_new+min(lambda)
   spectrum=INTERPOL_V8(spectrum_ln,lambda,lambda_linear,/NAN)
   error=INTERPOL_V8(error_ln,lambda,lambda_linear,/NAN)
   Mask=INTERPOL_V8(Mask_ln,lambda,lambda_linear,/NAN)

    new_Data=fltarr(NAXIS1_new,5)
	;oplotObj1 = Obj_New('cgOverPlot',lambda_linear,spectrum,thick=1, color="red") 
	;cgzplot,lambda,spectrum_ln,OPLOTS=oplotObj1

	new_Data[*,0]=spectrum
	new_Data[*,1]=error
	new_Data[*,2]=spectrum*0.0+!Values.F_NAN
	new_Data[*,3]=Mask ;mask from sdss
	new_Data[*,4]= 	spectrum*0.0+!Values.F_NAN  ; there is no sky spectrum in SDSS
	
	;Update the new header
		header=header_origin
		sxaddpar,header,'CRVAL1',min(lambda)
		sxaddpar,header,'NAXIS1',NAXIS1_new
		sxaddpar,header,'CD1_1',cdelt_new
		sxaddpar,header,'CRPIX1','0'
	    sxaddpar,header,'ROW1',"Spectrum"
	    sxaddpar,header,'ROW2',"Error"
	    sxaddpar,header,'ROW3',"Spectrum_nocalib"
	 	sxaddpar,header,'ROW4',"Error_nocalib"
	 	sxaddpar,header,'ROW5',"Sky"
 		sxaddpar,header,'SPECID',OBJECT_NAME
 		
	    mwrfits,new_Data,FILE_OUT,header,/Silent,/Create

return,1

END

PRO convert_GAMA, file, fileout
	 	
	Spec_struc=GAMA_READ_FITS2(file)
	Original_spec=Spec_struc.spectrum
	
	;Uniform the spectral resolution
	Spec_struc=GAMA_Uniform_resolution(Spec_struc)
	object=Spec_struc.name
	
	Sensibility_function=fltarr(n_elements(Spec_struc.Spectrum))
	Sensibility_function(*)=1.
	Spectrum_Corrected= Spec_struc.Spectrum 
	
	
 ;Write in fits file
 	DATA_final=fltarr(n_elements(Spec_struc.lambda),6)
 	DATA_final[*,0]=Spectrum_Corrected
 	Data_final[*,1]=Spec_struc.error
 	Data_final[*,2]=Spec_struc.Spectrum_nocalib
 	Data_final[*,3]=Spec_struc.Error_nocalib
 	Data_final[*,4]=Spec_struc.Sky
 	Data_final[*,5]=Sensibility_function
	header_final=Spec_struc.header
    sxaddpar,header_final,'ROW6',"Sensibility"
    sxaddpar,header_final,'NAXIS2','6'
    sxaddpar,header_final,'FLUCOR','T'
    sxaddpar,header_final,'DLam','1.6'
        
    mwrfits,DATA_final,fileout,header_final,/Silent,/Create
	;endif

END


PRO Convert_DATA

catalogue = '/Users/hypatia/Documents/GAMA/MZ_dwarf/Sample/File_upload.csv' 
PATH = '/Users/hypatia/Documents/GAMA/MZ_dwarf/Sample/GAMA_Original/'
PATH_OUT = '/Users/hypatia/Documents/GAMA/MZ_dwarf/Sample/GAMA_convert/'

readcol, catalogue, format='(a,a)',specid,url

for i=0, n_elements(specid) -1 do begin
survey=STRMID(URL(i), 40, 4)
print,survey

IF STRMATCH(survey,"gama") then begin
namefile=STRMID(URL(i), 59, 20)
remchar,namefile,'"'
print,namefile
convert_GAMA, PATH+namefile, PATH_OUT+specid(i)+".fits" 
ENDIF ELSE begin 
namefile=STRMID(URL(i), 45, 25)
remchar,namefile,'"'
print,namefile
res=CONVERT_SDSS_2_GAMA( PATH+namefile , specid(i), PATH_OUT+specid(i)+".fits" , /SILENT)
ENDELSE
endfor



END

