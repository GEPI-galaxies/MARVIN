PRO GAMA_Print_Emi_Catalogue

@GAMA_Main_Config

readcol,"GAMMA_Emission_lines_model.dat",format='(a,f,i)',lines_name, lines_ref,lines_group,/silent

restore,PATH_CATALOGUE+saved_lines

n_targets=n_elements(targets)
header='"Name"'
quantities=['Flux_','ErrFlux_','Sigma_','FluxLim_','SNR_','QC_','EW_']
quantities_type=['f12.5','f12.5', 'f10.5','f12.5','f7.2','i','f12.5']
type_print='(1(a,:, ", ")'
for m=0, n_elements(quantities)-1 do begin
	for n=0, n_elements(lines_name)-1 do begin 
	header= header +',"'+quantities[m]+lines_name[n]+'"'
	endfor
	type_print=type_print+','+strcompress(string(n_elements(lines_name)),/REMOVE_ALL)+'('+quantities_type(m)+',:, ", ")'
endfor
type_print=type_print+')'
print, type_print
STRREPLACE, saved_lines, '.sav', '.dat'

openw, U,PATH_CATALOGUE+saved_lines ,/GET_LUN
printf,U,format='(a)', header
for i=0, n_targets -1 do begin 
	 Name= targets[i]
	 Lines_array=LINE_ALL[*,i+1]
	 Continuum=Lines_array[*].Background 
	 Center_gauss=Lines_array[*].Center 
	 Sigma_gauss=Lines_array[*].Width
	 Flux=Lines_array[*].Flux
	 Err_Flux=Lines_array[*].Err_Flux_gaussian
	 FluxLim=Lines_array[*].Flux_limit
	 SN_detection=Lines_array[*].SNR  
	 QC_index=Lines_array[*].QC_index
	 EW=Lines_array[*].EW
	printf,U,format=type_print, name,Flux,Err_Flux,Sigma_gauss,FluxLim, SN_detection,QC_index,EW

	
endfor

FREE_LUN, U 
print, "Catalogue saved in :",saved_lines

END
