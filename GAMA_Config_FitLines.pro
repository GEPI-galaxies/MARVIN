;
;Configuration file for the emission and absorption lines fitting routine
;     GAMA_FIT_LINES
;     GAMA_MEASURE_INDEX

c=299792.458d                                       
pi=!pi
NOMINAL_RESOLUTION=5 ; in [A] GAMA Spectral resolution 

;Input spectra files 
preffixe="S_"
suffixe=".fits"
path="../Data_analyse/cb07wavelet/"

;Options for the autosave
AUTO_SAVE_P=5 				; Frequency of the auto-save
saved_file="Fit_lines.sav"  ; Save file for emission line
saved_file_index="Fit_indexlines.sav"  ; Save file
path_figure="../Data_analyse/cb07wavelet/" ; path to save the plot files

;OPTION for the emission line fitting 
default_Line_window=60.		; in [A], define the window from both side of the line(s) to consider when fitting
							;Should be large enough to cover both line and background to calculate the noise level
SIGMA_MIN=1.0				; in [A], Minimun possible value of the sigma for emission line
SIGMA_MAX=10.				; in [A], Maxinum possible value of the sigma for emission line
SIGMA_DEFAULT=2.8 			; in [A], Starting value for the sigma value of the gaussian model
WAVELENGTH_TOLERANCE=5. 	; in [A], Tolerance on the value of the gaussian center relative to the reference value
SN_detection_threshold=5.0    ; S/N Threshold for the detection of emission lines
MAX_ERROR_FLUX=0.85         ; Maxinum relative error on the flux. Reject line above this Threshold


;OPTION for index measurements
LICK_INDEX=1			; Measure the index at Lick resolution
LICK_RESOLUTION=10		; Spectral resolution of the Lick indices
FLAG_FIT_MODEL=0		; Measure index on modeled stellar continum	-> 1
FLAG_FIT_EMI=1			; Fit emission line before calculing index in e type index
BP_MAX_FRACTION=0.3		; Maxinum fraction of BP allowed in the line
SMOOTH_4PLOT=1          ; SMOOTH factor for the spectrum in plots
PRINT_RESULTS=0 		;Print results on screen
SAVE_CATALOGUE=1	    ; 1 to save the result in a catalogue in $Index_catalogue
PLOT_SAVE=0				; Save a plot of each index estimation
Index_catalogue="Gama_index_fit_catalogue.csv"