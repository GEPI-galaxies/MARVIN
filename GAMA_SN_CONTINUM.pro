	;-----------------------------------------------------------------------
	; FUNCTION GAMA_SN_continum, lambda, spectrum, region_SN
	;
	; DESCRIPTION : Calculate the SNR of a spectrum in a given lambda region
	;				Use sigma clipping to calculate a robust mean (sigma_clip=5)
	; 		INPUT : lambda and flux array
	;				region_SN 2-element vector defining the lambda region e.g.[3400:4000]	
	;		RETURN : SNR
	; 	
	;----------------------------------------------------------------------
 
	
	FUNCTION GAMA_SN_CONTINUM, lambda, spectrum, region_SN
	 region_flux= spectrum(WHERE(lambda GE region_SN[0] and lambda LE region_SN[1]))
	 RESISTANT_Mean, region_flux, 5, Mean_region, Sigma_region
	 Result = Moment(region_flux,/NAN, MDEV=mdev_region)
	 return, [Mean_region/mdev_region]
	END
