FUNCTION GAMA_RUN_STARLIGHT, objet, spectrum_str,lambda_vector=lambda_vector, Kinematics_vector=Kinematics_vector, FIX_KINEMATICS=FIX_KINEMATICS,sufixe_step

@GAMA_Main_Config	
	
	Dim_spectrum=n_elements(spectrum_str.Lambda)
	; Prepare input file for starlight
	file_cxt=objet+sufixe_step+".cxt"
	openw, U,PATH_FIT+file_cxt ,/GET_LUN
	 	for j=0,Dim_spectrum-1 do begin
			 printf,U,format='(f,e,e,i)',spectrum_str.Lambda(j),spectrum_str.Spectrum(j),spectrum_str.Error(j),spectrum_str.Mask(j)
		endfor
	FREE_LUN, U 	
	print, "        -Input file: "+file_cxt
	
	; Configure the grid for starlight

	IF  KEYWORD_SET(lambda_vector) then lambda_range=STRTRIM(lambda_vector,2) ELSE lambda_range=STRTRIM([min(spectrum_str.Lambda),max(spectrum_str.Lambda)],2)
	IF  KEYWORD_SET(Kinematics_vector) then kinematics=STRTRIM(Kinematics_vector,2) ELSE kinematics=['0','150']
	
	file_out_step=objet+sufixe_step+".BN"
	file_grid=objet+sufixe_step+".in"

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

	Comd_lambda1="'s/$Lambda_MIN_FIT/"+lambda_range(0)+"/' "
	Comd_lambda2="'s/$Lambda_MAX_FIT/"+lambda_range(1)+"/' "
	IF KEYWORD_SET(FIX_KINEMATICS) then Comd_fitkine="'s/$KIN_FIT/FXK/' " ELSE Comd_fitkine="'s/$KIN_FIT/FIT/' "
	Comd_input_file="'s/$INPUT/"+file_cxt+"/' "
	Comd_output_file="'s/$OUTPUT/"+file_out_step+"/' "
	Comd_vo="'s/$V0_STARLIGHT/"+Kinematics(0)+"/' "
	Comd_vd="'s/$VD_STARLIGHT/"+Kinematics(1)+"/' "
	
	command_line="sed -e "+Comd_config+"-e "+Comd_lambda1+"-e "+Comd_lambda2+"-e "+Comd_fitkine+"-e "+Comd_input_file+"-e "+Comd_output_file+"-e "+Comd_vo+"-e "+Comd_vd+ "< "+PATH_SOFT+"grid_file.in > "+ PATH_FIT+file_grid
	spawn,command_line
	print, "        -Input grid: "+file_grid
	
	; Run starlight	
	file_log_step=objet+sufixe_step+".log"
	start_starlight_step=PATH_STARLIGHT+"StarlightChains_v04.exe < "+PATH_FIT+file_grid	
	tofile= PATH_FIT+file_log_step
	spawn,start_starlight_step,listen
	
	result_fit_exist=file_exist(PATH_FIT+file_out_step)
	IF result_fit_exist then print,"        Starlight fitting: done"
	IF ~result_fit_exist then print,"        Starlight fitting: failed"
	return, result_fit_exist
		
END
