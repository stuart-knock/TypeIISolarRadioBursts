;Dynamic_Spectra.pro
;
; PURPOSE: This IDL script plots ripple parameters and dynamic 
;          spectra from the data output from Dynamic_Spectra.f90               
;
;
; INPUT: Inputs are read from files in outputDir which end with
;        fileident.
;
; OUTPUT: If plotting to file (ie selecting ps or eps at prompt) then 
;         then 3 files called Parameters1.ps, Parameters2.ps, and 
;         DS*.ps where *=Observer
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;--------------------------------------------------------------------
;***************** Enter Initialisation Stuff Here  *****************
;====================================================================

 fileident='temp.dat'
 outputDir='../OUTPUT/'  ;../OUTPUT/

;t_steps must be less then or equal to the number of time steps calculated
; <= wc -l Radial_Height_*  where *=fileident
 t_steps=173    ;Parameters: the number of time steps calculated
 Observer='3'  ;which observer; usually 1 | 2 | 3  
 Band='FH'     ;which emission bands to plot F|H|HF|FH (last 2 are equivalent)
;1sfu=10^-22 W m^-2 Hz^-1
 Background=[0,-22.0]; [0=no|1=yes, Log_10 lowest flux/background flux level]
 TIME_RANGE=[0,0,33 ]; [0=no|1=yes, UT_Begin(hours), UT_End(hours)]
 FREQ_RANGE=[0,4.0e4,5.0e8]; [0=no|1=yes, FREQ_low(Hz), FREQ_high(Hz)]
;********************************************************************
;______________END Enter Initialisation Stuff Here____________________


 xscale=1.0e9 
 yscale=1.0e9

 sol_rad=6.965e8    ;Solar radius in meters
 c=2.9979e8         ;speed of light
 AU=1.49599e11      ;1 Astronomical Unit
 Sr=6.965e8         ;1 Solar Radii
 Sp=2194560.0       ;1 Solar period
 Bo=1.89743128e-4   ;B-field strength at sun;  1AU by Parker model
 mu0=4.0*!pi /1.0e7
 mi=1.672614e-27    ;mass of ion (proton assumed)

;load 1AU Vsw  
 openr, lun, outputDir+'parameters_'+fileident, /get_lun
 readf, lun, FORMAT='(24(/),4X,E11.3)', Vsw_1AU
 free_lun,lun

;load array sizes
 openr, lun, outputDir+'LOG_'+fileident, /get_lun
 readf, lun, FORMAT='(9(/),30X,I4,2(/),40X,I4,2(/),32X,I5,2(/),28X,I4,12(/),36X,I4)', $ 
             Freq_Bins, time_bins_inw, time_inw, t_res, par_arr_size
 free_lun,lun

 t_res=t_res/3600.0 ;time resolution in hours 

;load how many time bins each time step corresponds to 
 openr, lun, outputDir+'scaling_'+fileident, /get_lun
 time=fltarr(t_steps)
 readf, lun, time
 free_lun,lun

 FOR j=1,t_steps-1 DO BEGIN
   time[j]=time[j-1]+time[j]
 ENDFOR
 time=time-time[0]
 time=time*t_res  ;time axis array for parameters
;print, time

 t_bins=LONG(time[t_steps-1]/t_res) ;Number of time bins in Dynamic Spectra 

;--------------------------------------------------------------------
 
 time_spectra=FIndGen(t_bins)*t_res ;time axis array for spectra

 IF(TIME_RANGE[0] EQ 0)THEN BEGIN 
   UT_Begin=0                    ;begin & 
   UT_End=time_spectra[t_bins-1] ;end times of calculations/plots in hours
   print, 'Full time range being used'
   print, UT_Begin,':',UT_End
 ENDIF ELSE BEGIN
   UT_Begin=TIME_RANGE[1] 
   UT_End=TIME_RANGE[2]  
   print, 'Entire time range of data not being used'
 ENDELSE
;--------------------------------------------------------------------

;load frequency ranges  
 openr, lun, outputDir+'parameters_'+fileident, /get_lun
 readf, lun, FORMAT='(26(/),8X,F4.1,/,10X,F4.1,/,6X,F4.1)', $ 
        exp_min, exp_range, l_val
 free_lun,lun

 Frequency=l_val*10.0^(exp_range*FIndGen(freq_bins)/(freq_bins-1)+exp_min)

 IF(FREQ_RANGE[0] EQ 0)THEN BEGIN 
   FREQ_low=l_val*10.0^(exp_min)
   FREQ_high=l_val*10.0^(exp_range+exp_min)
   print, 'Full frequency range being used'
 ENDIF ELSE BEGIN
   FREQ_low=FREQ_RANGE[1] 
   FREQ_high=FREQ_RANGE[2]  
   print, 'Entire frequency range of data not being used'
 ENDELSE
;--------------------------------------------------------------------

;read in all three observer locations 
 openr, lun, outputDir+'parameters_'+fileident, /get_lun
 readf, lun,  FORMAT= $ 
  '(32(/),2X,E11.2,/,2X,E12.2,3(/),2X,E11.2,/,2X,E12.2,3(/),2X,E11.2,/,2X,E12.2)', $
   Obs1X, Obs1Y, Obs2X, Obs2Y, Obs3X, Obs3Y
 free_lun,lun

 CASE Observer OF
   '1': begin
          ObserverX=LONG(Obs1X/xscale)
          ObserverY=LONG(Obs1Y/yscale)
        end
   '2': begin
          ObserverX=LONG(Obs2X/xscale)
          ObserverY=LONG(Obs2Y/yscale)
        end
   '3': begin
          ObserverX=LONG(Obs3X/xscale)
          ObserverY=LONG(Obs3Y/yscale)
        end
 ELSE: begin
         print, 'That Observer does not exist, sorry.'
         stop
       end
 ENDCASE
;print, ObserverX, ObserverY
;--------------------------------------------------------------------

;--------------------------------------------------------------------
;******************* Initialize output device 1**********************
;====================================================================

 opfilename='Parameters1'
 NHeaderRows=0
 x_cm=16.0 
 y_cm=23.0
 NCLev=7
 NGLev=8*NCLev

 SelectDevice, OPFilename, NHeaderRows, x_cm, y_cm, OPDevice, NCLev, NGLev
 IF(OPDevice ne 'q')THEN  begin
  
  ;multiple plots on page [0, across page, down page]
   !p.multi=[0,1,5]

  ;OVER-RIDE  plot margins
   !y.margin=[0,0]   ; [bottom,top], def=[4,2]
   !x.margin=[7,1]   ; [left,right]
  
  ;OVER-RIDE  object margins
   !y.omargin=[3,1]  ; [bottom,top], def=[4,2]
   !x.omargin=[0,5]  ; [left,right]
     
   CASE OPDevice OF
     'ps': begin
            TITLE_FONT='!4'
          end
     'eps': begin
            TITLE_FONT='!4'
          end
     'x': begin
            TITLE_FONT='!3'
          end
   ELSE: begin
           print, 'You have not set up the title font for this device'
           stop
         end
   ENDCASE 
  ;------------------------------------------------------------------
  

  ;--------------------------------------------------------------------
  ;***************************** Load Data ****************************
  ;--------------------------------------------------------------------

   openr, lun, outputDir+'UCME_'+fileident, /get_lun
   UCME=fltarr(t_steps)
   readf, lun, UCME
   free_lun,lun
  ;--------------------------------------------------------------------

   openr, lun, outputDir+'Radial_Height_'+fileident, /get_lun
   Radial_Height=fltarr(t_steps)
   readf, lun, Radial_Height
   free_lun,lun
  ;--------------------------------------------------------------------

   openr, lun, outputDir+'theta_'+fileident, /get_lun
   theta=fltarr(par_arr_size,t_steps)
   readf, lun, theta 
   free_lun,lun

   theta=ROTATE(theta,4)

  ;theta ranges for contours
   Maxtheta=Max(theta[WHERE(theta NE 0)],MIN=Mintheta)  
   theta_Range=Maxtheta-Mintheta
   print, 'Maxtheta= ',Maxtheta
   print, 'Mintheta= ',Mintheta
  ;--------------------------------------------------------------------

   openr, lun, outputDir+'b_'+fileident, /get_lun
   b=fltarr(par_arr_size,t_steps)
   readf, lun, b
   free_lun,lun

   b=ROTATE(b,4)

   ;ripple curvature constant ranges for contours
   Maxb=Max(b[WHERE(b NE 0)],MIN=Minb) 
   b_Range=Maxb-Minb
   print, 'Maxb= ',Maxb
   print, 'Minb= ',Minb
  ;--------------------------------------------------------------------

   openr, lun, outputDir+'U_'+fileident, /get_lun
   U=fltarr(par_arr_size,t_steps)
   readf, lun,U 
   free_lun,lun

   U=ROTATE(U,4)

  ;flow speed ranges for contours
   MaxU=Max(U[WHERE(U NE 0)],MIN=MinU)  
   U_Range=MaxU-MinU
   print, 'MaxU= ',MaxU
   print, 'MinU= ',MinU
  ;--------------------------------------------------------------------

  ;--------------------------------------------------------------------
  ;******************************** Plot ******************************
  ;====================================================================

   ;_________________
   ;
   ; Variable Trends 
   ;_________________


  ;Radial height ranges for Y-axis
   MINRadial_Height=MIN(Radial_Height/yscale,MAX=MAXRadial_Height)

   Plot, time, Radial_Height/yscale, XTICKNAME=[REPLICATE(' ', 8)], charsize=2.3,  $ 
   ytitle=TITLE_FONT+'Radial Height!3 (Gm)', xrange=[UT_Begin,UT_End], $
   yRange=[MINRadial_Height,MAXRadial_Height], $ ;, YLOG=1
   YTickFormat='(I3)'
  ;--------------------------------------------------------------------

   UCME=UCME/1.0e6

  ;CME/global shock speed ranges for Y-axis 
   MINUCME=MIN(UCME,MAX=MAXUCME)

   Plot, time, UCME, XTICKNAME=[REPLICATE(' ', 8)], charsize=2.3,  $ 
   ytitle=TITLE_FONT+'Speed of Shock!3 (Mm s!U-1!N)', xrange=[UT_Begin,UT_End], $
   yRange=[MINUCME,MAXUCME], $;, YLOG=1
   YTickFormat='(F3.1)'
  ;--------------------------------------------------------------------

  ;Load theta ranges for Y-axis
   openr, lun, outputDir+'LOG_'+fileident, /get_lun
   readf, lun, FORMAT='(35(/),9X,F13.8,/,9X,F13.8)', MAX_theta, MIN_theta
   free_lun,lun 

   MAX_theta=180.0D * MAX_theta / !pi
   MIN_theta=180.0D * MIN_theta / !pi

  ;Set greyscale/colourscale 
   stepg=theta_Range/(NGLev);
   grey =indgen(NGLev)*stepg+mintheta;

   contour, theta, time,                                                       $ 
     FINDGEN(par_arr_size)*(MAX_theta-MIN_theta)/par_arr_size + MIN_theta,     $ 
     YRANGE=[MIN_theta,MAX_theta], xrange=[UT_Begin,UT_End],                   $ 
     ytitle=TITLE_FONT+'!9q!3!DUB!N !3(!Uo!N)', XTICKNAME=[REPLICATE(' ', 8)], $ 
     levels=grey, /fill,c_colors=2+indgen(NGLev), charsize=2.3

   COLORBAR, POSITION = [0.43, 0.933, 0.59, 0.961],        $ 
             RANGE = [mintheta,maxtheta],                  $ 
             NCOLORS = NGLev,  BOTTOM = 3, CHARSIZE = 1.7, $ 
             FORMAT = '(I3)', DIVISIONS = 5, VERTICAL = 1, RIGHT = 1 
  ;--------------------------------------------------------------------

  ;Load ripple curvature constant ranges for Y-axis
   openr, lun, outputDir+'LOG_'+fileident, /get_lun
   readf, lun, FORMAT='(41(/),5X,F6.1,/,5X,F7.1)', MAX_b, MIN_b
   free_lun,lun

  ;Set greyscale/colourscale 
   stepg=b_Range/(NGLev);
   grey =indgen(NGLev)*stepg+minb;

   CONTOUR, b, time, FINDGEN(par_arr_size)*(MAX_b-MIN_b)/par_arr_size + MIN_b, $ 
            YRANGE=[MIN_b,MAX_b], xrange=[UT_Begin,UT_End],                    $ 
            ytitle=TITLE_FONT+'Ripple Curvature !3(Log[m!U-1!N])',             $ 
            XTICKNAME=[REPLICATE(' ', 8)],  charsize=2.3,                      $ 
            levels=grey, /fill,c_colors=2+indgen(NGLev)

   COLORBAR, POSITION = [0.24, 0.933, 0.40, 0.961],       $ 
             RANGE = [minb,maxb], NCOLORS = NGLev,        $ 
             BOTTOM = 3, CHARSIZE = 1.7, FORMAT = '(I3)', $ 
             DIVISIONS = 5, VERTICAL = 1, RIGHT = 1 
  ;--------------------------------------------------------------------

  ;Load flow speed ranges for Y-axis
   openr, lun, outputDir+'LOG_'+fileident, /get_lun
   readf, lun, FORMAT='(38(/),5X,E10.1,/,5X,E10.1)', MAX_U, MIN_U
   free_lun,lun

   ;Set greyscale/colourscale 
   stepg=U_Range/(NGLev);
   grey =indgen(NGLev)*stepg+minU;

   MAX_U=MAX_U/1.0e6
   MIN_U=MIN_U/1.0e6

   CONTOUR, U, time, FINDGEN(par_arr_size)*(MAX_U-MIN_U)/par_arr_size + MIN_U,  $ 
            YRANGE=[MIN_U,MAX_U], xrange=[UT_Begin,UT_End],                     $ 
            ytitle=TITLE_FONT+'Flow Speed !3(Mm s!U-1!N)',                      $ 
            xtitle=TITLE_FONT+'time!3 (hours)' ,                                $
            levels=grey, /fill,c_colors=2+indgen(NGLev), charsize=2.3

   COLORBAR, POSITION = [0.05, 0.933, 0.21, 0.961],       $ 
             RANGE = [minU,maxU], NCOLORS = NGLev,        $ 
             BOTTOM = 3, CHARSIZE = 1.7, FORMAT = '(I3)', $ 
             DIVISIONS = 5, VERTICAL = 1, RIGHT = 1 
  ;--------------------------------------------------------------------

 ENDIF ;(OPDevice ne 'q')
 CloseDevice, OPFilename, OPDevice
;__________________________________________________________________


;--------------------------------------------------------------------
;******************* Initialize output device 2**********************
;====================================================================

 opfilename='Parameters2'
 NHeaderRows=0
 x_cm=16.0 
 y_cm=23.0
 NCLev=7 
 NGLev=8*NCLev

 SelectDevice, OPFilename, NHeaderRows, x_cm, y_cm, OPDevice, NCLev, NGLev

 IF(OPDevice ne 'q')THEN  begin

   !p.multi=[0,1,4]

  ;OVER-RIDE  plot margins
   !y.margin=[0,0]   ; plot margins, [bottom,top], def=[4,2]
   !x.margin=[7,1]   ; [left,right]

  ;OVER-RIDE  object margins
   !y.omargin=[3,1]  ; object margins, [bottom,top], def=[4,2]
   !x.omargin=[0,5]  ; [left,right]


   CASE OPDevice OF
     'ps': begin
            TITLE_FONT='!4'
          end
     'eps': begin
            TITLE_FONT='!4'
          end
     'x': begin
            TITLE_FONT='!3'
          end
   ELSE: begin
           print, 'You have not set up the title font for this device'
           stop
         end
   ENDCASE 

  ;--------------------------------------------------------------------
  ;***************************** Load Data ****************************
  ;--------------------------------------------------------------------

   openr, lun, outputDir+'Ti_'+fileident, /get_lun
   Ti=fltarr(par_arr_size,t_steps)
   readf, lun, Ti
   free_lun,lun

   MaxTi=Max(Ti[WHERE(Ti NE 0)],MIN=MinTi)  
   Ti_Range=MaxTi-MinTi
   print, 'MaxTi= ',MaxTi
   print, 'MinTi= ',MinTi

   Ti=ROTATE(Ti,4)
  ;--------------------------------------------------------------------

   openr, lun, outputDir+'Te_'+fileident, /get_lun
   Te=fltarr(par_arr_size,t_steps)
   readf, lun, Te
   free_lun,lun

   MaxTe=Max(Te[WHERE(Te NE 0)],MIN=MinTe)  
   Te_Range=MaxTe-MinTe
   print, 'MaxTe= ',MaxTe
   print, 'MinTe= ',MinTe

   Te=ROTATE(Te,4)
  ;--------------------------------------------------------------------

   openr, lun, outputDir+'Ne_'+fileident, /get_lun
   N_e=fltarr(par_arr_size,t_steps)
   readf, lun, N_e
   free_lun,lun

   MaxN_e=Max(N_e[WHERE(N_e NE 0)],MIN=MinN_e)  
   N_e_Range=MaxN_e-MinN_e
   print, 'MaxN_e= ',MaxN_e
   print, 'MinN_e= ',MinN_e

   N_e=ROTATE(N_e,4)
  ;--------------------------------------------------------------------

   openr, lun, outputDir+'B1_'+fileident, /get_lun
   B1=fltarr(par_arr_size,t_steps)
   readf, lun, B1
   free_lun,lun

   MaxB1=Max(B1[WHERE(B1 NE 0)],MIN=MinB1)  
   B1_Range=MaxB1-MinB1
   print, 'MaxB1= ',MaxB1
   print, 'MinB1= ',MinB1

   B1=ROTATE(B1,4)
  ;--------------------------------------------------------------------


  ;--------------------------------------------------------------------
  ;******************************** Plot ******************************
  ;====================================================================

   ;_________________
   ;
   ; Variable Trends 
   ;_________________

   openr, lun, outputDir+'LOG_'+fileident, /get_lun
   readf, lun, FORMAT='(47(/),6X,F6.1,/,6X,F5.1)', MAX_N_e, MIN_N_e
   free_lun,lun

   print, Min_N_e, MAX_N_e

  ;Set greyscale/colourscale 
   stepg=N_e_Range/(NGLev);
   grey =indgen(NGLev)*stepg+minN_e;

   CONTOUR, N_e, time,                                                       $ 
            FINDGEN(par_arr_size)*(MAX_N_e-MIN_N_e)/par_arr_size + MIN_N_e,  $ 
            YRANGE=[MIN_N_e,MAX_N_e], xrange=[UT_Begin,UT_End],              $ 
            ytitle=TITLE_FONT+' N!De!N !3 (Log[m!U-3!N])',                   $ 
             XTICKNAME=[REPLICATE(' ', 8)],                                  $ 
            levels=grey, /fill,c_colors=2+indgen(NGLev), charsize=2.3

   COLORBAR, POSITION = [0.77, 0.933, 0.97, 0.961],       $ 
             RANGE = [minN_e,maxN_e], NCOLORS = NGLev,    $ 
             BOTTOM = 3, CHARSIZE = 1.7, FORMAT = '(I3)', $ 
             DIVISIONS = 5, VERTICAL = 1, RIGHT = 1 
  ;--------------------------------------------------------------------


   openr, lun, outputDir+'LOG_'+fileident, /get_lun
   readf, lun, FORMAT='(50(/),6X,E7.1,/,6X,E7.1)', MAX_Te, MIN_Te
   free_lun,lun

  ;Set greyscale/colourscale  
   stepg=Te_Range/(NGLev);
   grey =indgen(NGLev)*stepg+minTe;

   MAX_Te=MAX_Te/1.0e6
   MIN_Te=MIN_Te/1.0e6

   CONTOUR, Te, time, FINDGEN(par_arr_size)*(MAX_Te-MIN_Te)/par_arr_size + MIN_Te, $ 
            YRANGE=[MIN_Te,MAX_Te], xrange=[UT_Begin,UT_End],                      $ 
            ytitle=TITLE_FONT+' T!De!N !3(10!U6!N K)',                             $ 
            XTICKNAME=[REPLICATE(' ', 8)],                                         $ 
            levels=grey, /fill,c_colors=2+indgen(NGLev), charsize=2.3


   COLORBAR, POSITION = [0.535, 0.933, 0.735, 0.961],     $ 
             RANGE = [minTe,maxTe], NCOLORS = NGLev,      $ 
             BOTTOM = 3, CHARSIZE = 1.7, FORMAT = '(I3)', $ 
             DIVISIONS = 5, VERTICAL = 1, RIGHT = 1 
  ;--------------------------------------------------------------------


   openr, lun, outputDir+'LOG_'+fileident, /get_lun
   readf, lun, FORMAT='(53(/),6X,E7.1,/,6X,E7.1)', MAX_Ti, MIN_Ti
   free_lun,lun

  ;Set greyscale/colourscale 
   stepg=Ti_Range/(NGLev);
   grey =indgen(NGLev)*stepg+minTi;

   MAX_Ti=MAX_Ti/1.0e6
   MIN_Ti=MIN_Ti/1.0e6

   CONTOUR, Ti, time, FINDGEN(par_arr_size)*(MAX_Ti-MIN_Ti)/par_arr_size + MIN_Ti, $ 
            YRANGE=[MIN_Ti,MAX_Ti], xrange=[UT_Begin,UT_End],                      $ 
            ytitle=TITLE_FONT+' T!Di!N !3(10!U6!N K)',                             $ 
            XTICKNAME=[REPLICATE(' ', 8)],                                         $ 
            levels=grey, /fill,c_colors=2+indgen(NGLev), charsize=2.3


   COLORBAR, POSITION = [0.30, 0.933, 0.50, 0.961],       $ 
             RANGE = [minTi,maxTi], NCOLORS = NGLev,      $ 
             BOTTOM = 3, CHARSIZE = 1.7, FORMAT = '(I3)', $ 
             DIVISIONS = 5, VERTICAL = 1, RIGHT = 1 
  ;--------------------------------------------------------------------


   openr, lun, outputDir+'LOG_'+fileident, /get_lun
   readf, lun, FORMAT='(44(/),6X,F6.1,/,6X,F6.1)', MAX_B1, MIN_B1
   free_lun,lun

  ;Set greyscale/colourscale 
   stepg=B1_Range/(NGLev);
   grey =indgen(NGLev)*stepg+minB1;

   CONTOUR, B1, time, FINDGEN(par_arr_size)*(MAX_B1-MIN_B1)/par_arr_size + MIN_B1, $ 
            YRANGE=[MIN_B1,MAX_B1], xrange=[UT_Begin,UT_End],                      $ 
            ytitle=TITLE_FONT+' B!D1!N !3(Log[T])',                                $ 
            xtitle=TITLE_FONT+' time !3(hours)',                                   $
            levels=grey, /fill,c_colors=2+indgen(NGLev), charsize=2.3

   COLORBAR, POSITION = [0.065, 0.933, 0.265, 0.961],     $ 
             RANGE = [minB1,maxB1], NCOLORS = NGLev,      $ 
             BOTTOM = 3, CHARSIZE = 1.7, FORMAT = '(I3)', $ 
             DIVISIONS = 5, VERTICAL = 1, RIGHT = 1 
  ;--------------------------------------------------------------------

 ENDIF   ;(OPDevice ne 'q')
 CloseDevice, OPFilename, OPDevice

;__________________________________________________________________
;__________________________________________________________________


;--------------------------------------------------------------------
;******************* Initialize output device 3**********************
;====================================================================

 opfilename='DS'+Observer
 NHeaderRows=0
 x_cm=27.0 
 y_cm=17.0
 NCLev=7 
 NGLev=8*NCLev

 SelectDevice, OPFilename, NHeaderRows, x_cm, y_cm, OPDevice, NCLev, NGLev

 IF(OPDevice ne 'q')THEN  begin

   !p.multi=[0,1,1]

  ;OVER-RIDE  object margins
   !x.omargin=[0,0]     ; [left,right]
   !y.omargin=[0,0]     ; [bottom,top]


   CASE OPDevice OF
     'ps': begin
            TITLE_FONT='!4'
          end
     'eps': begin
            TITLE_FONT='!4'
          end
     'x': begin
            TITLE_FONT='!3'
          end
   ELSE: begin
           print, 'You have not set up the title font for this device'
           stop
         end
   ENDCASE 

  ;--------------------------------------------------------------------
  ;***************************** Load Data ****************************
  ;--------------------------------------------------------------------

   isitloaded=SIZE(Dynamic_SpectraF)
   print, isitloaded

   IF(SIZE(previousObs, /TYPE) EQ 0)THEN previousObs='0'
   IF(SIZE(previousrun, /TYPE) EQ 0)THEN previousrun='No previous run'

   IF((isitloaded[2] NE freq_bins) OR (isitloaded[1] NE  t_bins) $ 
      OR (previousObs NE Observer) OR (previousrun NE fileident))THEN BEGIN
     print, 'LOADING FUNDAMENTAL...'
     previousObs=Observer
     previousrun=fileident
     openr, lun, outputDir+'DS'+Observer+'F'+'_'+fileident, /get_lun
     Dynamic_SpectraF=fltarr(freq_bins,t_bins)
     readf, lun, Dynamic_SpectraF
     free_lun,lun

     minFluxF=MIN(Dynamic_SpectraF)
     print, 'minFluxF=', minFluxF

     Dynamic_SpectraF=ROTATE(Dynamic_SpectraF,4)
   ENDIF
  ;--------------------------------------------------------------------

   isitloaded=SIZE(Dynamic_SpectraH)
   print, isitloaded

   IF(SIZE(previousObs, /TYPE) EQ 0)THEN previousObs='0'
   IF(SIZE(previousrun, /TYPE) EQ 0)THEN previousrun='No previous run'

   IF((isitloaded[2] NE freq_bins) OR (isitloaded[1] NE  t_bins) $ 
      OR (previousObs NE Observer) OR (previousrun NE fileident))THEN BEGIN
     print, 'LOADING HARMONIC...'
     previousObs=Observer
     previousrun=fileident
     openr, lun, outputDir+'DS'+Observer+'H'+'_'+fileident, /get_lun
     Dynamic_SpectraH=fltarr(freq_bins,t_bins)
     readf, lun, Dynamic_SpectraH
     free_lun,lun

     minFluxH=MIN(Dynamic_SpectraH)
     print, 'minFluxH=', minFluxH

     Dynamic_SpectraH=ROTATE(Dynamic_SpectraH,4)
   ENDIF
  ;--------------------------------------------------------------------

  ;--------------------------------------------------------------------
  ;******************* Load colours & set ranges **********************
  ;====================================================================

  print, 'SETTING RANGES...'
  ;get stored colour table
   restore, 'colours/DS_W.sav'
   tvlct, r, g, b 

   CASE Band OF
     'H': begin
            Dynamic_Spectra=Dynamic_SpectraH 
            minFlux=minFluxH
          end
     'F': begin
            Dynamic_Spectra=Dynamic_SpectraF
            minFlux=minFluxF
          end
     'HF': begin
             DSF=10^Dynamic_SpectraF &  DSF=DSF-(DSF EQ 1)
             DSH=10^Dynamic_SpectraH &  DSH=DSH-(DSH EQ 1)
             Dynamic_Spectra=DSF+DSH
             Dynamic_Spectra[WHERE(Dynamic_Spectra NE 0.0)]= $ 
                  ALOG10(Dynamic_Spectra[WHERE(Dynamic_Spectra NE 0.0)])
             minFlux=MIN([minFluxF,minFluxH])
             print, minFlux
           end
     'FH': begin
             DSF=10^Dynamic_SpectraF &  DSF=DSF-(DSF EQ 1)
             DSH=10^Dynamic_SpectraH &  DSH=DSH-(DSH EQ 1)
             Dynamic_Spectra=DSF+DSH
             Dynamic_Spectra[WHERE(Dynamic_Spectra NE 0.0)]= $ 
                  ALOG10(Dynamic_Spectra[WHERE(Dynamic_Spectra NE 0.0)])
             minFlux=MIN([minFluxF,minFluxH])
             print, minFlux
           end
   ELSE: begin
           print, 'Band should be one of [F,H,HF,FH]'
           stop
         end
   ENDCASE 

   Dynamic_Spectra=(Dynamic_Spectra EQ 0.0)*minFlux + TEMPORARY(Dynamic_Spectra)

  ;Arbitrarily define lower Flux plot value 
   IF(Background[0] EQ 1)THEN BEGIN 
     D_Spectra=Dynamic_Spectra > (Background[1]) 
     minFlux = Background[1]
     print, 'minFlux reset to ', minFlux
   ENDIF ELSE BEGIN 
     D_Spectra=Dynamic_Spectra
   ENDELSE
   maxFlux=MAX(D_Spectra, MIN=minFlux)
   print, 'maxFlux=', maxFlux 

  ;Set greyscale/colourscale 
   stepg=(maxFlux-minFlux)/(NGLev-1);
   grey =indgen(NGLev)*stepg+minFlux;

  ;--------------------------------------------------------------------
  ;******************************** Plot ******************************
  ;====================================================================


  ;OVER-RIDE  plot margins

   !x.margin=[7,1]     ; [left,right]
   !y.margin=[8,2]     ; [bottom,top]

   ;________________
   ;
   ; Dynamic Spectra
   ;________________

  print, 'PLOTTING...' 

   contour, D_Spectra, time_spectra, Frequency, $ 
   title=TITLE_FONT+'Dynamic Spectra !3 Observer at ('+STRTRIM(STRING(ObserverX),2)+','+STRTRIM(STRING(ObserverY),2)+') Gm', $
   xtitle=TITLE_FONT+'time !3 (hours)', ytitle=TITLE_FONT+'Frequency !3 (Hz)',  $ 
   YLOG=1,  xrange=[UT_Begin,UT_End], yrange=[FREQ_low,FREQ_high], $ 
   levels=grey, /fill,c_colors=1+indgen(NGLev), charsize=1.3
  ;--------------------------------------------------------------------
  ;Frequency drift rates

   Average_Ne=(1.36e6*(Radial_Height/Sol_rad)^(-2.14D) $    ;SAITO POLAND MUNRO 1977
              +1.68e8*(Radial_Height/Sol_rad)^(-6.13D))*1.0e6

   Average_fp=8.976D*SQRT(Average_Ne)    
   ObserverX=ObserverX*xscale    &   ObserverY=ObserverY*yscale
   OPLOT, time+(SQRT((ObserverX-Radial_Height)^2+(ObserverY)^2)/c)/3600 , $ 
          Average_fp,  THICK=2 
   OPLOT, time+(SQRT((ObserverX-Radial_Height)^2+(ObserverY)^2)/c)/3600 , $ 
          2*Average_fp,  THICK=2 
  ;--------------------------------------------------------------------       

   ;______________
   ;
   ; Add colourbar 
   ;______________ 

   COLORBAR, TITLE = 'log!D10!N '+TITLE_FONT+'Flux !3(W m!U-2!N Hz!U-1!N Sr!U-1!N)',  $ 
   POSITION = [0.07, 0.03, 0.93, 0.06],  NCOLORS = NGLev-1, $ 
   RANGE = [minFlux, maxFlux], BOTTOM = 2, CHARSIZE = 1.1,  $ 
   FORMAT = '(F6.2)', DIVISIONS = 5

;__________________________________________________________________
 ENDIF   ;(OPDevice ne 'q')
 CloseDevice, OPFilename, OPDevice
;__________________________________________________________________
;__________________________________________________________________

end ;Dynamic_Spectra.pro  
