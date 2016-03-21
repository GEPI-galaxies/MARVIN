	;-----------------------------------------------------------------------
	; PROCEDURE GAMA_ZOOM_PLOT,Lambda,spectrum, Plot_definition,oplotObj=oplotObj
	;
	; DESCRIPTION : Zoom into a window cgplot using the mouse cursor
	; 				Need the cgplot library
	; 		INPUT  : Lambda and spectrum array
	;				 Object with the array of vectors to overplot
	;				 Plot_definition structure
	;					-> Plot_definition{ Title, ytitle, charsize, xtitle}
	;----------------------------------------------------------------------
	
	PRO GAMA_ZOOM_PLOT, Lambda,spectrum, Plot_definition,oplotObj=oplotObj
				cursor,x1,y1,/DATA,/DOWN
				cursor,x2,y2,/DATA,/UP
				x_zoom=[min([x1,x2]),max([x1,x2])]
				y_zoom=[min([y1,y2]),max([y1,y2])]			
	
				IF KEYWORD_SET(oplotObj) THEN BEGIN 
				cgPlot, Lambda,spectrum, OPLOTS=oplotObj $
					,xrange=x_zoom,yrange=y_zoom,XS=1,Title=Plot_definition.Title,ytitle=Plot_definition.ytitle,charsi=Plot_definition.charsize,$
						xtitle=Plot_definition.xtitle
				ENDIF ELSE begin
				cgPlot, Lambda,spectrum,xrange=x_zoom,yrange=y_zoom,XS=1,Title=Plot_definition.Title,ytitle=Plot_definition.ytitle,$
				charsi=Plot_definition.charsize,xtitle=Plot_definition.xtitle
				ENDELSE
	END

