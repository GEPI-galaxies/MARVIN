;-----------------------------------------------------------------------
	; FUNCTION Read_GAMA_spectrum
	;
	; DESCRIPTION : LOAD GAMA fit files in a structure 
	; 				Read raw GAMA files and GAMA_soft processed files
	; 		INPUT  : file path
	;		OUTPUT : structur
	;				 	lambda:lambda array
	;					spectrum:spectrum array
	;					error:error array
	;					sky:sky array
	;					name:OBJECT_NAME from header
	; 					header: header from the file
	;----------------------------------------------------------------------
	
	FUNCTION GAMA_READ_FITS,file,SILENT=SILENT
	
		if (FILE_EXIST(file)) then begin
		Data=readfits(file,header_origin, EXTEN_NO=0, /SILENT)
	 	;Check the data size 
	 	check_data=size(DATA)
			IF (check_data(0) NE 2) then begin 
			 	IF ~(KEYWORD_SET(SILENT)) then print,"Wrong fits format"
			return,0
			endif	
		ENDIF ELSE begin 
				IF ~(KEYWORD_SET(SILENT)) then print,"File not found"
		return,0
		ENDELSE
	
	    OBJECT_NAME=SXPAR(header_origin,'SPECID')
  		CRVAL=SXPAR(header_origin,'CRVAL1')
	 	NAXIS1=SXPAR(header_origin,'NAXIS1')
	 	CDELT=SXPAR(header_origin,'CD1_1')
	 	CRPIX=SXPAR(header_origin,'CRPIX1')
	 	wstart = CRVAL + (1 - CRPIX) * CDELT
	 	lambda=indgen(NAXIS1)*CDELT+wstart


		Data_tags = SXPAR( header_origin, 'ROW*',COUNT=n_rows)
		Data_tags=STRCOMPRESS(Data_tags,/REMOVE_ALL)
		Data_tags=['Lambda',Data_tags,'name']
		tag_descript="D("+STRCOMPRESS(NAXIS1,/REMOVE_ALL)+")"
		for i=0, n_rows -1 do begin 
		tag_descript=tag_descript+","+"D("+STRCOMPRESS(NAXIS1,/REMOVE_ALL)+")"
		endfor
		tag_descript=tag_descript+",A(1)"
		CREATE_STRUCT, spectrum_struc,'' ,Data_tags,  tag_descript
		ADD_TAG,spectrum_struc, 'Header',header_origin ,  Final_struc
		Final_struc.(0)=lambda
		Final_struc.name=OBJECT_NAME
		for i=1, n_rows  do begin
		Final_struc.(i)=Data(*,i-1)
		endfor

	 	return,Final_struc

	END
