;Calculate the mean sensibility correction to apply to the spectrum of a field

PRO GAMA_Sensibility_Field, Configuration

@GAMA_Main_Config

List_Spec = FILE_SEARCH(PATH_SENS+'Fit/','Sens_'+Configuration+'*', count=n_file)

Lambda_resample=indgen(8860-3725) + 3725
Dim_spectrum=n_elements(Lambda_resample)
	
Sensibility_total= fltarr(Dim_spectrum, n_elements(List_Spec))
position_fiber=fltarr(2,n_elements(List_Spec))

cgPS_Open, PATH_SENS+'/Fig/'+Configuration+'_sens.ps'
cgplot,[-1,-1],[-1,1], title=configuration+" - "+TRIM(string(n_elements(List_Spec)),2)+' fibers',xrange=[3740,8800.],XS=1,yrange=[0,2],xtitle='Wavelength (Angstrom)',ytitle='Sensivity function [spectroscopy/models]', charthi=2, thick=2

; LoadCT, 33, /Silent
;color_vector=BytScl(indgen(n_elements(List_Spec)))

For i=0, n_file -1 do begin
print,List_Spec(i)
	Spec_struc=GAMA_READ_FITS(List_spec(i))
	check_data=size(Spec_struc)
	IF (check_data(0) EQ 0) then begin 
		print, "        !! Failed to open file "+List_spec(i)		
	ENDIF ELSE begin 
        Sensibility_total[*,i]=INTERPOL_V8(Spec_struc.Sensibility,Spec_struc.Lambda,Lambda_resample)
	cgoplot,Spec_struc.Lambda,Spec_struc.Sensibility,Color='Grey'
	ENDELSE	
				

	
endfor 


	Sigma_Sens_configuration=fltarr(Dim_spectrum)
	Sens_configuration=fltarr(Dim_spectrum)
	Quantil1=fltarr(Dim_spectrum)
	Quantil2=fltarr(Dim_spectrum)
	for j=0, Dim_spectrum -1 do begin
	summary=Quartils(Sensibility_total[j,*])
	Sens_configuration[j]=summary[2]
	Quantil1[j]=summary[1]
	Quantil2[j]=summary[3]
	endfor
	cgoplot,Lambda_resample,replicate(1,n_elements(Lambda_resample)),color='Charcoal',thick=1 
	cgoplot,Lambda_resample,Gauss_smooth(Sens_configuration,300, /EDGE_MIRROR),color='red',thick=2
	cgoplot,Lambda_resample,Quantil1,color='red',thick=2, linesty=2
	cgoplot,Lambda_resample,Quantil2,color='red',thick=2, linesty=2 
cgPS_Close, /PNG

DATA_final=fltarr(Dim_spectrum,3)
DATA_final[*,0]=Gauss_smooth(Sens_configuration,300, /EDGE_MIRROR)
DATA_final[*,1]=Quantil1
DATA_final[*,2]=Quantil2

MKHDR, header_sens, DATA_final
sxaddpar,header_sens,'CRVAL1','3725'
sxaddpar,header_sens,'NAXIS1',Dim_spectrum
sxaddpar,header_sens,'CD1_1','1'
sxaddpar,header_sens,'CRPIX1','0'
WRITEFITS,PATH_SENS+'/Sens_function/'+configuration+'_Sens.fits',DATA_final,header

END


PRO GAMA_Sensibility_Field_ALL
@GAMA_Main_Config

list_file= FILE_SEARCH(PATH_SENS+'/Fit/','Sens*',count=n_files)
Config=strarr(n_elements(list_file))

for i=0d, n_files -1 do begin
Config[i]=STRMID(list_file[i], STRLEN(PATH_SENS+'/Fit/')+4, 10)
endfor
Config = Config[ rem_dup(Config)]

for i=0d, n_elements(Config) -1 do begin
GAMA_Sensibility_Field, Config[i]
endfor
END
