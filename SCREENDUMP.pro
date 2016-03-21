  PRO SCREENDUMP, filename, Encapsulated=encapsulated

      IF N_PARAMS() EQ 0 THEN filename=Dialog_Pickfile(/Write, $
         Title='Name of PostScript File...', File='screendump.ps')

      ; Get the screen dump of the current graphics window.
      screenDump = cgSnapShot()

      ; Make a PostScript window with the same aspect
      ; ratio as the current display window. Use color
      ; PostScript if the COLOR keyword is set.
      aspect = PSWINDOW()

      ; Open a PostScript file and dump it.
      thisDevice = !D.NAME
      SET_PLOT, 'PS', /Copy
      DEVICE, FILENAME=filename, XSIZE=aspect.xsize, YSIZE=aspect.ysize, $
         XOFFSET=aspect.xoffset, YOFFSET=aspect.yoffset, COLOR=1, $
         ENCAPSULATED=Keyword_Set(encapsulated), Inches=aspect.inches, $
         DECOMPOSED=1

      ; Display the screen dump. Fill up the window.
      cgImage, screenDump
      DEVICE, /CLOSE_FILE
      SET_PLOT, thisDevice

   END