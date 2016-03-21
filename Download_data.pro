PRO Download_data

catalogue = 'File_upload.csv' 

PATH_LOCAL_DATA = '/Users/hypatia/Documents/GAMA/Data/'
PATH_OUT = '/Users/hypatia/Documents/GAMA/MZ_dwarf/Sample/GAMA_Original/'

readcol, catalogue, format='(a,a)',specid,url
openw, lun, 'download_missing.sh', /get_lu

for i=0, n_elements(specid) -1 do begin
namefile=STRMID(URL(i), 60, 20)
remchar,namefile,'"'

IF (FILE_EXIST(PATH_LOCAL_DATA+namefile)) then begin
SPAWN,'cp '+PATH_LOCAL_DATA+namefile+' '+PATH_OUT+namefile
ENDIF ELSE begin 
printf,lun,format='(a)', 'curl -O --user gama:aaomega '+URL(i)
ENDELSE
endfor

 free_lun, lun

END