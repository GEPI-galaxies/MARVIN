

PRO GAMA_QF_stat, Q_Flag
print," --------- Statistic ----------"
print,""
class=WHERE(STRCMP(Q_Flag,'X'),n_class)
print, "[Unclassified]= "+STRTRIM(n_class,2)
class=WHERE(STRCMP(Q_Flag,'L'),n_class)
print, "[Pending]= "+STRTRIM(n_class,2)
print,""
class=WHERE(STRCMP(Q_Flag,'G'),n_class)
print, "[Quality ok]= "+STRTRIM(n_class,2)
print,""
class=WHERE(STRCMP(Q_Flag,'F'),n_class)
print, "[Frinding]= "+STRTRIM(n_class,2)
class=WHERE(STRCMP(Q_Flag,'W'),n_class)
print, "[Wrong detection]= "+STRTRIM(n_class,2)
class=WHERE(STRCMP(Q_Flag,'S'),n_class)
print, "[Bad splicing]= "+STRTRIM(n_class,2)
class=WHERE(STRCMP(Q_Flag,'C'),n_class)
print, "[Bad flux calibration]= "+STRTRIM(n_class,2)
class=WHERE(STRCMP(Q_Flag,'A'),n_class)
print, "[AGN]= "+STRTRIM(n_class,2)
print,""
print, "[TOTAL]= "+STRTRIM(n_elements(Q_Flag),2)
print, "------------------------------"
END

PRO GAMA_QF_Help
print," --------- Help for GAMA Quality routine ----------"
	print, " "
	print,"         * Interative options *  "
	print,'z -> Zoom (Select area to zoom)'	
	print,'r --> redraw the plot'
	print,'p --> go to previous object'
	print,'space --> print statistic'
	print, ""
	print, "         * Quality flags *        "
	print,'n --> Good quality and go to next target'
	print,'w --> wrong detection'
	print,	'f --> frindging'
	print,	's --> bad splicing'
	print,	'c --> bad flux calibration'
	print,	'a --> AGN with broad emission'
	print,	'l --> classify later on'
	print, "------------------------------"

END 

;-------------------------------------------------------------
; Routine to create an interactive window for quality control	
;
;	INPUT : Lambda, Spectrum and additional arrays of the spectrum to be plot
;			Plot_definition-> Structure with the option of the plot
;
; 	OUTPUT : QF_obj -> Quality flag given to the spectrum 
;-------------------------------------------------------------

PRO GAMA_Quality_interactive, Lambda,spectrum,oplotObj,Plot_definition, QF_obj	
	

	cgPlot, Lambda,spectrum, OPLOTS=oplotObj, Title=Plot_definition.Title, $
			ytitle=Plot_definition.ytitle,xtitle=Plot_definition.xtitle, $
			charsi=Plot_definition.charsize,yrange=Plot_definition.yrange,YS=1,xrange=minmax(Lambda),xs=1
	
	interactive_loop=1
		while interactive_loop EQ 1 do begin
		action_event=GET_KBRD(1) 
		if STRCMP(action_event,'z') then GAMA_ZOOM_PLOT,Lambda,spectrum,Plot_definition,oplotObj=oplotObj
		if STRCMP(action_event,'h') then GAMA_QF_Help ; Print help
	    if STRCMP(action_event,'p') then begin 
									QF_obj='P' ;Go to previous
									interactive_loop=0
								    endif		
		if STRCMP(action_event,' ') then begin 
									QF_obj=' ' ;PRINT statistic
									interactive_loop=0
								    endif		
		if STRCMP(action_event,'b') then begin 
									QF_obj='B' ;Show gama webpage
									interactive_loop=0
								    endif
		 	    					  
		if STRCMP(action_event,'n') then begin
									QF_obj='G' ;Good and go to next
									interactive_loop=0
								    endif
		if STRCMP(action_event,'r') then GAMA_REDRAW, Lambda,spectrum,Plot_definition,oplotObj=oplotObj
		if STRCMP(action_event,'w') then begin 
										QF_obj='W' ;wrong detection
										interactive_loop=0
								    endif
		if STRCMP(action_event,'f') then begin
									QF_obj='F' ;frindging
									interactive_loop=0
								    endif
		if STRCMP(action_event,'s') then begin
									QF_obj='S' ;bad splicing
									interactive_loop=0
								    endif
		if STRCMP(action_event,'c') then begin 
									QF_obj='C' ;bad flux calibration
									interactive_loop=0
								    endif
		if STRCMP(action_event,'a') then begin 
									QF_obj='A' ;AGN with broad emission
									interactive_loop=0
								    endif
		if STRCMP(action_event,'l') then begin 
									QF_obj='L' ; classify later on
									interactive_loop=0
								    endif
	    if STRCMP(action_event,'q') then begin 
									QF_obj='Q' ;Go to previous
									interactive_loop=0
								    endif
		endwhile
END	


;-----------------------------------------------------------------------
; ROUTINE GAMA_Quality_check
;
; DESCRIPTION : Routine to interactively perform the quality check of spectra
; 
;       INPUT : ascii file of a Catalogue with the format SPECID,CATAID,redshift 
;	 			The file is define in GAMA_Config_Quality.pro 
;	 			
;      OUTPUT : .sav file with the same format as the input plus a flag parameter
;     
;  INTERACTIVE OPTION : 
;
;  	 	z --> Zoom (Select area to zoom)	
;  		r --> redraw the plot
;  	 	p --> go to previous object
;		b --> show gama webpage of the object
;  	 	space --> print statistic
;  		 	* Quality flags *        
;		n --> Good quality and go to next target
;		w --> wrong detection
;		f --> frindging
;		s --> bad splicing
;		c --> bad flux calibration
;		a --> AGN with broad emission
;		l --> classify later on
;----------------------------------------------------------------------


PRO GAMA_Quality_check
@GAMA_Main_Config

	;Load catalogue defined in the config file
	catalogue_exist=file_exist(PATH_CATALOGUE+CATALOGUE_INPUT+catalogue_extension)
	case catalogue_exist of 
		'1': begin
			readcol, PATH_CATALOGUE+CATALOGUE_INPUT+catalogue_extension, format='(a,l,f,f,f,f,f,f)', SPECID,CATAID,redshift_gama, MAG_PETRO_U, MAG_PETRO_G, MAG_PETRO_R,  MAG_PETRO_I, MAG_PETRO_Z
			n_targets=n_elements(CATAID)
	 		end
		'0': begin
			 print, "Catologue "+CATALOGUE_INPUT+catalogue_extension+" didn't found in "+PATH_CATALOGUE
			 return
	end
	 
	endcase


; Main option menu
; Start a new classification (0) or load a saved classification
 
	print, "--------------------------------------------- "
	print, "-    GAMA spectrum quality control routine  -"
	print, "  "
	print, "Load Catalogue : ", CATALOGUE_INPUT+catalogue_extension
	print, "                with ", STRTRIM(n_targets,2), " targets"
	print, " "
	print, "--------------------------------------------- "
	print," What do you want to do?"
	print,"-> 0 Start a new classification"
	print,"-> 1 Load a saved classification"

loop_main_menu:
	action_event=GET_KBRD(1) 

	case action_event of 
  		'0': begin
			 print,''	
		  	 print," Starting classification from scratch"
		  	 print,"    "
		  	 n_end=n_targets
		  	 array_2_class=indgen(n_targets)
  			 Q_Flag=strarr(n_targets)
		  	 Q_Flag(*)="X"
		  	 SN_cont=fltarr(n_targets)
		  	 SN_cont(*)=-999.

		     saved_file=PREFIXE_QC0+CATALOGUE_INPUT+".sav"
		  		if file_exist(PATH_CATALOGUE+saved_file) then begin
		  			print, "File already exist, create a new file"
		       		saved_file=DIALOG_PICKFILE(/Write, FILTER = '*.sav',DEFAULT_EXTENSION='*.sav', PATH=PATH_CATALOGUE)
		        endif
		 	 save, SPECID,CATAID,redshift_gama, Q_Flag,SN_cont, filename=PATH_CATALOGUE+saved_file
			 print,"The classification will be save in "+PATH_CATALOGUE+saved_file
		  end

  		'1': begin
   			 print,"Load classification from file"
   			 saved_file=DIALOG_PICKFILE(/Read, FILTER = '*.sav', PATH=PATH_CATALOGUE)
			 restore,saved_file
			 tmp=WHERE(STRCMP(Q_Flag,'X') OR STRCMP(Q_Flag,'M') ,noclass)
			 tmp=WHERE(STRCMP(Q_Flag,'L'),Pending)
			 loop_saved_menu:
   		   	 print,"--> 0  Classify all unclassified target 'X' ("+STRTRIM(noclass,2)+' obj)'
   			 print,"--> 1  Classify pending target 'L' ("+STRTRIM(Pending,2)+' obj)'
			 option_load=GET_KBRD(1)
				case option_load of
					'0': begin 
					 array_2_class=WHERE(STRCMP(Q_Flag,'X') OR STRCMP(Q_Flag,'M'),n_end) 
					 end
					'1': begin 
					 array_2_class=WHERE(STRCMP(Q_Flag,'L'),n_end) 
					 end
					 else: begin
			 		   print,'Say again. Valid option are 0/1' 
			 		   goto, loop_saved_menu 
			 		 end
				endcase
	  		 end	  
   			     
	 	else:   begin 
			print,'Say again' 
			goto, loop_main_menu
			end
		endcase


; Start the interactive classification 
	print, ""
	print," Go to the interactive window .."

for n=0L, n_end-1 do begin 

	i=array_2_class(n)
	redshift=redshift_gama(i)
	spectrum_file=SPECID(i)+".fits"
print,PATH_OBS+spectrum_file
	Data_struct=GAMA_READ_FITS(PATH_OBS+spectrum_file)	
	check_data=size(Data_struct)
	
	IF ~(check_data(0) EQ 0) then begin
	spectrum=Data_struct.Spectrum
	smoothed_spec=smooth(spectrum, 6, /NAN)
	lambda=Data_struct.lambda
	error=Data_struct.error
	sky=Data_struct.sky
	
	;SN_Broad=SXPAR(header_origin,'SN_ALT') 		; Median S/N per pixel in 4500-6500 A band   
	;SN_runz=SXPAR(header_origin,'SN') 				; S/N as measured by runz 
	SN=[GAMA_SN_continum(lambda/(redshift+1), spectrum, BLUE_REGION),GAMA_SN_continum(lambda, spectrum, RED_REGION)]

	;Spectro-photometry
	mag_photo=[MAG_PETRO_U(i),MAG_PETRO_G(i),MAG_PETRO_R(i),MAG_PETRO_I(i),MAG_PETRO_Z(i)]

	info_spectrophoto=GAMA_specphoto(lambda,spectrum,mag_photo)


;Plot options
	bound=minmax(smoothed_spec(WHERE(lambda GT min(lambda)+50 and lambda LT max(lambda)-50)),/NAN)
	Yscale=[bound[0]*0.8,bound[1]*1.2]
	Plot_definition={Title:SPECID[i]+" SNR="+STRTRIM(SN[0],2)+'/'+STRTRIM(SN[1],2)+" Z="$
					+STRTRIM(redshift_gama(i),2), xtitle:"Lambda [A]",ytitle:'Flux 1E-17 ergs/s/cm2/A)',charsize:1.2, yrange:Yscale}
	oplotObj1= Obj_New('cgOverPlot',Lambda,sky/10., color='BLU5')
	oplotObj2= Obj_New('cgOverPlot',[5700,5700],[-1000,1000], thick=2, color='red') ; Dicrohic region
	oplotObj3= Obj_New('cgOverPlot',info_spectrophoto.lambda_filters, info_spectrophoto.flux_photo_filters, psym=6,symsize=2, color='red') ; Spectrophotometry_catalogue
	oplotObj4= Obj_New('cgOverPlot',info_spectrophoto.lambda_filters, info_spectrophoto.flux_spec, psym=16,symsize=2, color='grey') ; Spectrophotometry_spectrum

	oplotObj=[oplotObj1,oplotObj2,oplotObj3,oplotObj4]
	
	object_lines=[3727.,4340.,4862., 4959.,5006.8,6300.,6562.,5755.,6717.,6731.]*(1+redshift_gama(i))
			for l=0, n_elements(object_lines)-1 do begin 
			tmp_obj= Obj_New('cgOverPlot',REPLICATE(object_lines(l),2),[-1000.,1000.],linestyle=2, thick=1,color='TG6')
			oplotObj=[oplotObj,tmp_obj]
			endfor
	object_lines=[3933.,3968.,5889.]*(1+redshift_gama(i))
			for l=0, n_elements(object_lines)-1 do begin 
			tmp_obj= Obj_New('cgOverPlot',REPLICATE(object_lines(l),2),[-1000.,1000.],linestyle=2, thick=1,color='PUR4')
			oplotObj=[oplotObj,tmp_obj]
			endfor

	GAMA_Quality_interactive, Lambda,smoothed_spec,oplotObj,Plot_definition, QF_obj
		case QF_obj of
		'P': begin
			n=n-2
			Q_Flag(i)='X'
			end
		'Q': begin
			Q_Flag(i)='X'
			save, SPECID,CATAID,redshift_gama, Q_Flag,SN_cont, filename=saved_file
			print,"Classification interrupted at "+SPECID(i)
			print,"Saved in:", saved_file
			GAMA_QF_stat, Q_Flag
			return
			end
		' ': begin	
			n=n-1
			GAMA_QF_stat, Q_Flag
		    end
	    'B': begin
			Q_Flag(i)='X'
			open_url=url_Path+SPECID(i)+url_suffixe
			mg_open_url, open_url
			n=n-1
			end

	 	else: begin Q_Flag(i)=QF_obj
			print,Q_Flag(i)
			end		
		endcase
	
		ENDIF ELSE Q_Flag(i)="M" ; missing file
			;Auto-save
			IF ((i MOD AUTO_SAVE_QC0) EQ 0) then begin
			print,SPECID(i),CATAID(i), Q_Flag(i)
			save, SPECID,CATAID,redshift_gama, Q_Flag,SN_cont, filename=saved_file
			print,"Saved in:", saved_file
			ENDIF

endfor
save, SPECID,CATAID,redshift_gama, Q_Flag,SN_cont, filename=saved_file
print,"!! Congratulation !! Classification finished"
print,"Saved in:", saved_file
GAMA_QF_stat, Q_Flag

END



