;-----------------------------------------------------------------------
; ROUTINE PRINT_QC_Catalogue
;
; DESCRIPTION : Save the quality check catalogue in .sav format to .dat
;
;----------------------------------------------------------------------

PRO PRINT_QC0_Catalogue
@GAMA_Main_Config

saved_file=DIALOG_PICKFILE(/Write, FILTER = '*.sav',DEFAULT_EXTENSION='*.sav', PATH=PATH_CATALOGUE)

restore, saved_file
n_targets=n_elements(CATAID)
STRREPLACE,saved_file,".sav",".dat"

openw, U,saved_file ,/GET_LUN
for i=0L, n_targets-1 do begin
print,Q_FLAG(i)
printf,U,format='(a,a,i,a,f,a,a)', SPECID(i)," , ",CATAID(i)," , ",redshift_gama(i)," , ",Q_FLAG(i)
endfor
FREE_LUN, U 
print, "Catalogue saved in :",saved_file
END
