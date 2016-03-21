
FUNCTION GAMA_SN_RN, lambda,spectrum, fit, error_stat,mask

	residual_fit= spectrum-fit

	SN_window_down=[4500,5400,6000,6800]
	SN_window_up=[4700,5600, 6200,7000]

	rN_sN_window=fltarr(4)
	rN_sN_total=-999
	rN_sN_window=[!Values.F_NAN,!Values.F_NAN,!Values.F_NAN,!Values.F_NAN]
	S_sN_window=fltarr(4)
	S_sN_window=[!Values.F_NAN,!Values.F_NAN,!Values.F_NAN,!Values.F_NAN]
	S_sN_total=-999

		FOR kk=0, 3 do begin
		region=WHERE(lambda GE SN_window_down(kk) and lambda LE SN_window_up(kk) and mask EQ 0, n_region)
		IF n_region GT 1 then begin
		level_error=mean(error_stat(region),/NAN)
		level_residual=sigma(residual_fit(region))

		 rN_sN=level_residual/level_error
		 S_sN=mean(spectrum(region),/NAN)/level_error
		 
		 rN_sN_window(kk)=rN_sN 
		 S_sN_window(kk)=S_sN
		ENDIF
		ENDFOR

	arN_sn=max(rN_sN_window,tmp,/NAN)
	REMOVE,tmp ,rN_sN_window
	aS_sn=max(S_sN_window,tmp2,/NAN)
	REMOVE, tmp2,S_sN_window

	rN_sN_total= MEAN(rN_sN_window,/NAN);+0.3
	S_sN_total= MEAN(S_sN_window,/NAN)
	return, [rN_sN_total,S_sN_total]

END
