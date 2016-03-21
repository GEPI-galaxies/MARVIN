
;-------------------------------------------------------------
; Routine to convert the saved quality control .sav into .dat	
;-----------------------------------------------------------
PRO GAMA_QC1_Catalogue
@GAMA_Main_Config

restore, filename=QC1_CATALOGUE_OUTPUT+".sav"
 
n_targets=n_elements(cataid)

openw, U,QC1_CATALOGUE_OUTPUT+".dat" ,/GET_LUN
printf,U,format='(a)'," SPECID  CATAID  z  QC_Flag    QF Q1_Flag rN_sN  S_sN"
for i=0L, n_targets-1 do begin
printf,U,format='(a,2x,i,2x,f8.6,2x,a,2x,i,2x,a,2x,f7.4,2x,f7.4)', spectrumid(i),cataid(i),z(i),QC_Flag(i),QF(i),Q1_Flag(i),rN_sN(i),S_sN(i)
endfor
FREE_LUN, U 
print, "Catalogue saved in :",QC1_CATALOGUE_OUTPUT+".dat" 

END





;-------------------------------------------------------------
; Routine to create an interactive window for quality control	
;
;	INPUT : Lambda, Spectrum and additional arrays of the spectrum to be plot
;			Plot_definition-> Structure with the option of the plot
;
; 	OUTPUT : QF_obj -> Quality flag given to the spectrum 
;-------------------------------------------------------------

PRO GAMA_Fit_check, Lambda,spectrum,oplotObj,Plot_definition, QF_obj	
	

	cgPlot, Lambda,spectrum, OPLOTS=oplotObj, Title=Plot_definition.Title, $
			ytitle=Plot_definition.ytitle,xtitle=Plot_definition.xtitle, $
			charsi=Plot_definition.charsize,yrange=Plot_definition.yrange,YS=1,xrange=minmax(Lambda),xs=1
	
	interactive_loop=1
		while interactive_loop EQ 1 do begin
		action_event=GET_KBRD(1) 
		if STRCMP(action_event,'z') then GAMA_ZOOM_PLOT,Lambda,spectrum,Plot_definition,oplotObj=oplotObj
		if STRCMP(action_event,'p') then begin 
									QF_obj='P' ;Go to previous
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
; ROUTINE GAMA_Fit_check_all
;
; DESCRIPTION : Routine to interactively perform the quality check of stellar continmum fit 
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
;  		 	* Quality flags *        
;		n --> Good quality and go to next target
;		w --> wrong detection
;		l --> classify later on
;----------------------------------------------------------------------


PRO GAMA_Fit_check_all

@GAMA_Main_Config

	;Load catalogue defined in the config file
	catalogue_exist=file_exist(QC1_CATALOGUE_INPUT)
	case catalogue_exist of 
		'1': begin
			restore,QC1_CATALOGUE_INPUT
			n_targets=n_elements(CATAID)
	 		end
		'0': begin
			 print, "Catologue "+QC1_CATALOGUE_INPUT+" didn't found in "+PATH_FIT
			 return
			end
	 
	endcase


	


; Main option menu
; Start a new classification (0) or load a saved classification
 
	print, "--------------------------------------------- "
	print, "-    GAMA Starlight fit quality control routine  -"
	print, "  "
	print, "Load Catalogue : ", QC1_CATALOGUE_INPUT
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
			 ;Prepare QF1
			 Q1_Flag=strarr(n_targets)
			 Q1_Flag(*)="-9"
			 tmp=WHERE(QF EQ 1 OR QF EQ -4 ,TBD)
			 Q1_Flag(tmp)="X"

		  	 rN_sN=fltarr(n_targets)
		  	 rN_sN(*)=!Values.F_NAN
 			 S_sN=fltarr(n_targets)
		  	 S_sN(*)=!Values.F_NAN

		     	output_sav=QC1_CATALOGUE_OUTPUT+'.sav'
print,output_sav
		  		if file_exist(output_sav) then begin
		  			print, "File already exist, create a new file"
		       		output_sav=DIALOG_PICKFILE(/Write, FILTER = '*.sav',DEFAULT_EXTENSION='*.sav', PATH=PATH_CATALOGUE)
		      	        endif
		 	 save, spectrumid,cataid,z, QC_Flag,QF,Q1_Flag,rN_sN,S_sN, filename=output_sav
			 print,"The classification will be save in "+output_sav
		 	 end

  		'1': begin
   			 print,"Load classification from file"
   			 output_sav=DIALOG_PICKFILE(/Read, FILTER = '*.sav', PATH=PATH_CATALOGUE)
			 restore,output_sav
			 array_2_class=WHERE(STRCMP(Q1_Flag,'X') OR STRCMP(Q1_Flag,'L') ,n_end)

   		   	 print,"--> Check all unclassified fit 'X' or 'L' ("+STRTRIM(n_end,2)+' obj)'
   			 
			end			
		

	else:   begin 
			print,'Say again' 
			goto, loop_main_menu
		end
	endcase

; Start the interactive classification 
	print, ""
	print," Go to the interactive window .."

FOR n=0L, n_end-1 do begin 

	i=array_2_class(n)
	spectrum_file=spectrumid(i)+".fits"
	Data_struct=GAMA_READ_FITS(PATH_FIT+preffixe_final+spectrum_file)
	check_data=size(Data_struct)
	
	IF (check_data(0) EQ 0) then begin
	Q1_Flag(i)="NF" 
	
	ENDIF ELSE begin	

	spectrum=Data_struct.Spectrum
	smoothed_spec=smooth(spectrum, 6, /NAN)
	lambda=Data_struct.lambda
	error=Data_struct.error
	mask=Data_struct.mask
	sky=Data_struct.sky
	Continuum=Data_struct.Continuum
	
	SN_RN=GAMA_SN_RN( lambda,spectrum, Continuum, error,mask)
	S_Sn(i)=SN_RN[1]
	rN_Sn(i)=SN_RN[0]
	SN=[GAMA_SN_continum(lambda, spectrum, BLUE_REGION),GAMA_SN_continum(lambda, spectrum, RED_REGION)]

;Plot options
	bound=minmax(smoothed_spec(WHERE(lambda GT min(lambda)+50 and lambda LT max(lambda)-50)),/NAN)
	Yscale=[bound[0]*0.8,bound[1]*1.2]
	Plot_definition={Title:spectrumid[i]+" SNR="+STRTRIM(SN[0],2)+'/'+STRTRIM(SN[1],2)+" Z="$
					+STRTRIM(z(i),2), xtitle:"Lambda [A]",ytitle:'Flux 1E-17 ergs/s/cm2/A)',charsize:1.2, yrange:Yscale}
	oplotObj1= Obj_New('cgOverPlot',Lambda,sky/50., color='BLU5')
	oplotObj2= Obj_New('cgOverPlot',Lambda,smooth(Continuum,6), thick=2, color='red') ; Dicrohic region
	oplotObj=[oplotObj1,oplotObj2]
	
	object_lines=[3727.,4340.,4862., 4959.,5006.8,6300.,6562.,5755.,6717.,6731.]
			for l=0, n_elements(object_lines)-1 do begin 
			tmp_obj= Obj_New('cgOverPlot',REPLICATE(object_lines(l),2),[-1000.,1000.],linestyle=2, thick=1,color='TG6')
			oplotObj=[oplotObj,tmp_obj]
			endfor
	object_lines=[3933.,3968.,5889.]
			for l=0, n_elements(object_lines)-1 do begin 
			tmp_obj= Obj_New('cgOverPlot',REPLICATE(object_lines(l),2),[-1000.,1000.],linestyle=2, thick=1,color='PUR4')
			oplotObj=[oplotObj,tmp_obj]
			endfor

	GAMA_Fit_check, Lambda,smoothed_spec,oplotObj,Plot_definition, QF_obj
		case QF_obj of
		'P': begin
			n=n-2
			Q1_Flag(i)='X'
			end
		'Q': begin
			Q1_Flag(i)='X'
			save,  spectrumid,cataid,z, QC_Flag,QF,Q1_Flag,rN_sN,S_sN , filename=output_sav
			print,"Classification interrupted at "+spectrumid(i)
			return
			end
		
	 
	 	else: begin 
			Q1_Flag(i)=QF_obj
			end
		endcase
	
	ENDELSE
	
       ;Auto-save	
	IF ((i MOD AUTO_SAVE_QC1) EQ 0) then begin
	save,  spectrumid,cataid,z, QC_Flag,QF,Q1_Flag,rN_sN,S_sN,filename=output_sav
	ENDIF


ENDFOR

print,"!! Congratulation !! Classification finished"

END

