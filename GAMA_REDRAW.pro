	;-----------------------------------------------------------------------
	; PROCEDURE GAMA_REDRAW,Lambda,spectrum, Plot_definition,oplotObj=oplotObj
	;
	; DESCRIPTION : Redraw into a window cgplot. E.g: Use after GAMA_Zoom
	; 				Need the cgplot library
	; 		INPUT  : Lambda and spectrum array
	;				 Object with the array of vectors to overplot
	;				 Plot_definition structure
	;					-> Plot_definition{ Title, ytitle, charsize, xtitle}
	;----------------------------------------------------------------------

	PRO GAMA_REDRAW, Lambda,spectrum,Plot_definition,oplotObj=oplotObj
		IF KEYWORD_SET(oplotObj) THEN BEGIN 
		cgPlot, Lambda,spectrum, OPLOTS=oplotObj,Title=Plot_definition.Title,ytitle=Plot_definition.ytitle,$
				charsi=Plot_definition.charsize,xtitle=Plot_definition.xtitle,$
				YS=1,xrange=minmax(Lambda),xs=1,yrange=Plot_definition.yrange
	ENDIF ELSE begin
		cgPlot, Lambda,spectrum,Title=Plot_definition.Title,ytitle=Plot_definition.ytitle,charsi=Plot_definition.charsize,$
				xtitle=Plot_definition.xtitle,YS=1,xrange=minmax(Lambda),xs=1,yrange=Plot_definition.yrange
	ENDELSE
	END

