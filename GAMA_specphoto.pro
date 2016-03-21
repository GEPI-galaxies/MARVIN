
FUNCTION GAMA_specphoto,lambda,spectrum,mag_photo

@GAMA_Main_Config


 ; SDSS filters u,g,r,i,z
 lambda_central=[3557.,4702.,6165.,7491.,8946.] ; from http://www.astro.ljmu.ac.uk/~ikb/research/mags-fluxes/
 lambda_width=[550.,1300.,1200.,1300.,1000.]
 lambda_ub=4118.50
 lambda_iz=8206.


 flux_photo= mag2flux(mag_photo, ABwave=lambda_central)

 no_photo=WHERE(mag_photo GE 99.,cnt_photo)
	if cnt_photo GE 1 then  flux_photo(no_photo)=!Values.F_NAN	
 flux_ub=mean([flux_photo[0],flux_photo[1]])
 flux_iz=mean([flux_photo[3],flux_photo[4]])
	
 lambda_filters=[lambda_central[0],lambda_ub,lambda_central[1],lambda_central[2],lambda_central[3],lambda_iz,lambda_central[4]]
 flux_photo_filters=[flux_photo[0],flux_ub,flux_photo[1],flux_photo[2],flux_photo[3],flux_iz,flux_photo[4]]/factor_normalization


 ;Smooth spectrum
 Spec_sliding=GAMA_SLIDING_MOMENTS(spectrum,60)
 smoothed_spec=Spec_sliding(*,0)
 Spec_U=MEAN(Spec_sliding(WHERE(lambda GT lambda_central[0]- (min(lambda)-20) and lambda LT lambda_central[0]+ lambda_width[1]/2. )))
 Spec_ub=MEDIAN(smoothed_spec(WHERE(lambda GT lambda_ub -10 and lambda LT lambda_ub +10 )))
 Spec_G=MEAN(Spec_sliding(WHERE(lambda GT lambda_central[1]- lambda_width[1]/2. and lambda LT lambda_central[1]+ lambda_width[1]/2. )))
 Spec_R=MEAN(Spec_sliding(WHERE(lambda GT lambda_central[2]- lambda_width[2]/2. and lambda LT lambda_central[2]+ lambda_width[2]/2. )))
 Spec_I=MEAN(Spec_sliding(WHERE(lambda GT lambda_central[3]- lambda_width[3]/2. and lambda LT lambda_central[3]+ lambda_width[3]/2. )))
 Spec_iz=MEDIAN(smoothed_spec(WHERE(lambda GT lambda_iz -10 and lambda LT lambda_iz +10 )))

 flux_spec=[Spec_U,Spec_ub,Spec_G,Spec_R,Spec_I,Spec_IZ]	
 ;Calculate the aperture correction using g-,r- and i- band
 tmp=flux_photo_filters/flux_spec
 Aperture_correction=median([tmp(2),tmp(3),tmp(4)]) ;calculated in R-band

 
;Diff in flux
 Dif_UB=(flux_spec[1]-flux_photo_filters[1]/Aperture_correction)/(flux_photo_filters[1]/Aperture_correction)
 Dif_G=(flux_spec[2]-flux_photo_filters[2]/Aperture_correction)/(flux_photo_filters[2]/Aperture_correction)
 Dif_R=(flux_spec[3]-flux_photo_filters[3]/Aperture_correction)/(flux_photo_filters[3]/Aperture_correction)
 Dif_I=(flux_spec[4]-flux_photo_filters[4]/Aperture_correction)/(flux_photo_filters[4]/Aperture_correction)
 Dif_IZ=(flux_spec[5]-flux_photo_filters[5]/Aperture_correction)/(flux_photo_filters[5]/Aperture_correction)
 ErrAll_spectrophoto=[Dif_UB,Dif_G,Dif_R,Dif_R,Dif_IZ]
 Err_spectrophoto=sqrt(Dif_G^2+ Dif_R^2+Dif_R^2)
print,"Spectrophotometry uncertainties[%]:"
print,"UB-band   G-band   R-band  IZ-band"
print,ErrAll_spectrophoto*100.
print, "[total:]"
print,Err_spectrophoto*100.

spectrophoto_info={flux_spec : flux_spec, flux_photo_filters: flux_photo_filters/Aperture_correction,lambda_filters:lambda_filters,Aperture_correction:Aperture_correction}

return, spectrophoto_info
END
