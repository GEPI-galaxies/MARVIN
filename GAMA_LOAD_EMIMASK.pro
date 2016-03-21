;-----------------------------------------------------------------------
	; PROCEDURE GAMA_LOAD_emiMask
	; DESCRIPTION : Add the parameter from the configuration file in the 
	;               starlight grid
	;---------------------------------------------------------------------
	
	FUNCTION GAMA_LOAD_EMIMASK, lambda
	
	@GAMA_Main_Config

	Mask_lines=lambda*0
	readcol, PATH_MASK+MASK_STARLIGHT, format='(f,f,f,a)', l_min, l_max, weight, comments, /SILENT
	for i=0, n_elements(l_min)-1 do begin
	region_mask=WHERE(lambda GE l_min(i) and lambda LT l_max(i), n_pixels)
	IF n_pixels GT 1 then Mask_lines(region_mask)=1.
	endfor
	return, Mask_lines
	END

