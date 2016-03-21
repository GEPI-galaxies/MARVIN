PRO GAMA_Verify_StarlightFit, object,rN_sN_total,S_sN_total, PLOT_FIT=PLOT_FIT, PLOT_FIGURE=PLOT_FIGURE
	
preffixe="S_"
suffixe=".fits"
path="/Users/marodrig/Documents/GAMA/Test/Out_starlight/"
file=preffixe+object+suffixe
Data=readfits(path+file,header_origin, EXTEN_NO=0, /SILENT)

check_data=size(DATA)
IF (check_data(0) EQ 2) then begin

spectrum=Data[*,0]
error=Data[*,1]
Sensibility=Data[*,2]
fit_normalized=Data[*,3]
sky=Data[*,4]
mask=Data[*,5]

CRVAL=SXPAR(header_origin,'CRVAL1')
NAXIS1=SXPAR(header_origin,'NAXIS1')
CDELT=SXPAR(header_origin,'CD1_1')
CRPIX=SXPAR(header_origin,'CRPIX1')
wstart = CRVAL + (1 - CRPIX) * CDELT
lambda=indgen(NAXIS1)*CDELT+wstart

SN_WINDOW=SXPAR(header_origin,'SN')

error_stat=(error-(spectrum*abs(sensibility-1)) ) ; subtract the error from the sensibility correction/Users/marodrig/Documents/GAMA/GAMA_soft/GAMA_Main_Config.pro
residual_fit=(spectrum-fit_normalized)
un_correct_spec=spectrum*Sensibility

IF KEYWORD_SET(PLOT_FIT) then begin
oplotObj1 = Obj_New('cgOverPlot',lambda,fit_normalized,thick=2, color="red") 
oplotObj2 = Obj_New('cgOverPlot',lambda,error_stat,thick=1, color="green") 
oplotObj3 = Obj_New('cgOverPlot',lambda,residual_fit,thick=1, color="violet") 
oplotobjs=[oplotObj1,oplotObj2,oplotObj3]
CGZPLOT, lambda,spectrum, OPLOT=oplotobjs, title=object
ENDIF

IF KEYWORD_SET(PLOT_FIGURE) then begin
PS_Start,path+'Stellar_continuum.png',/quiet

oplotObj1 = Obj_New('cgOverPlot',lambda,smooth(spectrum,3),thick=1, color="black") 
oplotObj2 = Obj_New('cgOverPlot',lambda,smooth(fit_normalized,3),thick=2, color="red") 
;oplotObj3 = Obj_New('cgOverPlot',lambda,error_stat,thick=1, color="green") 
oplotobjs=[oplotObj1,oplotObj2]
!X.OMargin = [0.05, 0.05]
   !Y.OMargin = [0.05, 0.05]
!P.Multi = [0, 1, 2]
CGPLOT, lambda,smooth(un_correct_spec,3), color='grey',thick=2, OPLOT=oplotobjs, title=object,xra=[3600,7700],yrange=[5,40],xs=1,ys=1,/YLOG, XTIT='Wavelength [A]', YTITLE='Normalized Flux'
oplotObj1 = Obj_New('cgOverPlot',lambda,lambda*0.,thick=1, color="blue") 

CGPLOT, lambda,smooth(residual_fit,3), color='black', OPLOT=oplotObj1,thick=1,xra=[3600,7700],yrange=[-5,30],xs=1,ys=1, XTIT='Wavelength [A]', YTITLE='Normalized Flux'
PS_End
ENDIF

SN_window_down=[4500,5400,6000,6800]
SN_window_up=[4700,5600, 6200,7000]

rN_sN_window=fltarr(4)
rN_sN_total=-999
rN_sN_window=[!Values.F_NAN,!Values.F_NAN,!Values.F_NAN,!Values.F_NAN]
S_sN_window=fltarr(4)
S_sN_window=[!Values.F_NAN,!Values.F_NAN,!Values.F_NAN,!Values.F_NAN]
S_sN_total=-999

for i=0, 3 do begin
region=WHERE(lambda GE SN_window_down(i) and lambda LE SN_window_up(i) and mask EQ 0, n_region)
 IF n_region GT 1 then begin
 level_error=mean(error_stat(region),/NAN)
 level_residual=sigma(residual_fit(region))

 rN_sN=level_residual/level_error
 S_sN=mean(spectrum(region),/NAN)/level_error
 
 rN_sN_window(i)=rN_sN 
 S_sN_window(i)=S_sN
endif
endfor
arN_sn=max(rN_sN_window,tmp,/NAN)
REMOVE,tmp ,rN_sN_window
aS_sn=max(S_sN_window,tmp2,/NAN)
REMOVE, tmp2,S_sN_window

rN_sN_total= MEAN(rN_sN_window,/NAN);+0.3
S_sN_total= MEAN(S_sN_window,/NAN)
ENDIF

END


PRO GAMA_Verify_Fit_all
path="/Users/marodrig/Documents/GAMA/Test/Out_starlight/"

readcol, path+"list",format='(a)',object
n_obj=n_elements(object)

rN_sN=fltarr(n_obj)
SN=fltarr(n_obj)
for i=0L, n_obj -1 do begin
target=object(i)
GAMA_Verify_StarlightFit,target ,rN_sN_i,SN_i
rN_sN(i)=rN_sN_i;+0.5
SN(i)=SN_i
endfor



xrange = minmax(SN)
 yrange = minmax(rN_sN)
   xbinsize = .5
   ybinsize = 0.20
   
   
density = Hist_2D(SN,rN_sN, Min1=xrange[0], Max1=xrange[1], Bin1=xbinsize, $
                           Min2=yrange[0], Max2=yrange[1], Bin2=ybinsize)   
                           
   maxDensity = Ceil(Max(density)/1e2) * 1e2
   scaledDensity = BytScl(density, Min=0, Max=maxDensity)

;PS_Start,path+'Quality_continum.pdf',/quiet

 ;cgDisplay                          
   ; Load the color table for the display. All zero values will be gray.
   cgLoadCT, 33
   TVLCT, cgColor('white', /Triple), 0
   TVLCT, r, g, b, /Get
   palette = [ [r], [g], [b] ]
   
   ; Display the density plot.

   cgImage, scaledDensity, XRange=xrange, YRange=yrange, /Axes, Palette=palette, $
      XTitle='S/sN', YTitle='rN/sN', $
      Position=[0.125, 0.125, 0.9, 0.8]

       cgColorbar, Position=[0.125, 0.875, 0.9, 0.915], Title='N# galaxies', $
       Range=[0, maxDensity], NColors=254, Bottom=1, OOB_Low='white', $
       TLocation='Top'
;PS_End
save, rN_sN, SN, object, filename='FitQuality_testdata.save'

END

PRO GAMA_Sensitivity_all

path="/Users/marodrig/Documents/GAMA/Test/Out_starlight/"

readcol, path+"list",format='(a)',object
n_obj=n_elements(object)
preffixe="S_"
suffixe=".fits"

N_spectral=5200
lambda_final=indgen(N_spectral)*1.0+3800.

Sensitivity_all=MAKE_ARRAY(N_spectral,n_obj,/float)
Sensitivity_all(*)=!Values.F_NAN
amount_correction=fltarr(n_obj)
amount_correction(*)=!Values.F_NAN


;FOR i=0,500 do begin
FOR i=0, n_obj-1 do begin

file=preffixe+object[i]+suffixe

Data=readfits(path+file,header_origin, EXTEN_NO=0, /SILENT)
redshift=SXPAR(header_origin,'Z')
CRVAL=SXPAR(header_origin,'CRVAL1')
NAXIS1=SXPAR(header_origin,'NAXIS1')
CDELT=SXPAR(header_origin,'CD1_1')
CRPIX=SXPAR(header_origin,'CRPIX1')
wstart = CRVAL + (1 - CRPIX) * CDELT
lambda_obs=(indgen(NAXIS1)*CDELT+wstart)*(1+redshift)

check_data=size(DATA)
	IF (check_data(0) EQ 2) then begin	
	spectrum=Data[*,0]
	error=Data[*,1]
	Sensibility=Data[*,2]
	Sensitivity_all[*,i]=INTERPOL(Sensibility,lambda_obs,lambda_final)
	amount_correction[i]=TOTAL(ABS(1-Sensitivity_all[*,i]),/NAN)/N_spectral
	ENDIF
ENDFOR


MEANABSDEV_sensibility=fltarr(N_spectral)
mean_sensibility=fltarr(N_spectral)
MEANABSDEV_sensibility(*)=!Values.F_NAN
mean_sensibility(*)=!Values.F_NAN

for i=0, N_spectral-1 do begin

;sigma_sensibility[i]=MEANABSDEV(Sensitivity_all[i,*],/NAN )
MEANABSDEV_sensibility[i]=MEANABSDEV(Sensitivity_all[i,*],/NAN )
mean_sensibility[i]=MEAN(Sensitivity_all[i,*],/NAN )

endfor
PS_Start,path+'Sensitivity_correction_simulated_data.ps',/quiet
cgplot,lambda_final,mean_sensibility,xtit="Lambda [A]",ytit="Sensitivity corrective function",yrange=[0.5,1.5],xrange=[3800,9000],ys=1,xs=1


cgPolygon,[lambda_final,REVERSE(lambda_final)],[mean_sensibility+1*MEANABSDEV_sensibility,REVERSE(mean_sensibility-1*MEANABSDEV_sensibility)], fcolor='BLK3' ,/fill, color='BLK6'
cgOplot,lambda_final,mean_sensibility
;cgOplot, lambda_final,mean_sensibility+1*MEANABSDEV_sensibility,linesty=2, color='bleu'
cgOplot, lambda_final,lambda_final*0.+1., color='blue', linestyle=3

PS_End
;cgplot,lambda_final,MEANABSDEV_sensibility/mean_sensibility
END