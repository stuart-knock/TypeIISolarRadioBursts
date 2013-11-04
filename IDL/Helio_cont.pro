;Helio_cont.pro  
;
; PURPOSE: This IDL script plots ripple parameters as a function 
;          of location from the data output from Dynamic_Spectra.f90 
;
; INPUT: Inputs are read from files in outputDir which end with
;        fileident.
;
; OUTPUT: If plotting to file (ie selecting ps or eps at prompt) then 
;         then 1 file called Helio*.ps where *=Observer is written.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;--------------------------------------------------------------------
;***************** Enter Initialisation Stuff Here  *****************
;====================================================================
 fileident='temp.dat'
 outputDir='../OUTPUT/'  ;../OUTPUT/

;t_steps can be determined from parameters.dat, the number of times UT is written
; <= wc -l Radial_Height_*  where *=fileident
 t_steps=173   ;Parameters: the number of time steps calculated
 Observer='3'  ;which observer; usually 1 | 2 | 3 
 Plot_RANGE=[0,1.0] ;Overide spatial range [0=no|1=yes, PlotBox (Gm)]
;********************************************************************
;______________END Enter Initialisation Stuff Here____________________

 xscale=1.0e9 
 yscale=1.0e9
 XH=33   ;
 YH=33   ; 
 k=100   ;plot every kth global shock
  
 temp=LONARR(XH,2*YH-1)

;constants
 pi=3.141592653589793238462643D 
 sol_rad=6.965e8   ;Solar radius in meters
 arr=100   ;
;--------------------------------------------------------------------

;Load array of radial heights at each time step:
 openr, hi, outputDir+'Radial_Height_'+fileident, /get_lun
 Radial_Height=Dblarr(t_steps)
 readf, hi, Radial_Height
 free_lun,hi

 IF(Plot_RANGE[0] EQ 0)THEN BEGIN 
   RangeGm=Radial_Height[t_steps-1]/1.0e9
 ENDIF ELSE BEGIN
   RangeGm=Plot_RANGE[1]  
   print, 'Entire range of data not being used'
 ENDELSE
;--------------------------------------------------------------------

;Load array of the number of sources for each time step:
 openr, num, outputDir+'Source_Num_'+fileident, /get_lun
 Number_of_Sources=fltarr(t_steps)
 readf, num, Number_of_Sources
 free_lun, num

 print, 'Total(Number_of_Sources)=',LONG(Total(Number_of_Sources))

;Load array Wind parameters as a function of location for each time step, 
;in the form x, y, theta, Ne, B1, Ti, Te, Vsw, Kappa
 openr, loc, outputDir+'Wind_Loc_'+fileident, /get_lun
 Sources = Dblarr(9, LONG(Total(Number_of_Sources)))
 readf, loc, Sources 
 free_lun, loc    

;Load array of the source location and strength for each time step, 
;in the form x, y, Fund_Flux, Harm_Flux:l
 openr, loc, outputDir+'Source_Loc'+Observer+'_'+fileident, /get_lun
 SourcesFH = Dblarr(5, LONG(Total(Number_of_Sources)))
 readf, loc, SourcesFH 
 free_lun, loc    

;Load array of radial heights at each time step:
 openr, hi, outputDir+'CME_b_'+fileident, /get_lun
 CMEb=fltarr(t_steps)
 readf, hi, CMEb
 free_lun,hi
;--------------------------------------------------------------------
   
;Obtain MAX & MIN nonzero [theta, N_e, B1, Ti, Te, Vsw, Kappa]: 
 Maxtheta=MAX(Sources[2,WHERE(Sources[2,*] NE 0.0)],MIN=Mintheta) 
 MaxNe=MAX(Sources[3,WHERE(Sources[3,*] NE 0.0)],MIN=MinNe) 
 MaxB1=MAX(Sources[4,WHERE(Sources[4,*] NE 0.0)],MIN=MinB1) 
 MaxTi=MAX(Sources[5,WHERE(Sources[5,*] NE 0.0)],MIN=MinTi)   
 MaxTe=MAX(Sources[6,WHERE(Sources[6,*] NE 0.0)],MIN=MinTe) 
 MaxVsw=MAX(Sources[7,WHERE(Sources[7,*] NE 0.0)],MIN=MinVsw) 
 Maxkappa=MAX(Sources[8,WHERE(Sources[8,*] NE 0.0)],MIN=Minkappa) 
 print, 'Mintheta=',Mintheta, '   Maxtheta=',Maxtheta
 print, 'MinNe=',MinNe, '   MaxNe=',MaxNe
 print, 'MinB1=',MinB1, '   MaxB1=',MaxB1
 print, 'MinTi=',MinTi, '   MaxTi=',MaxTi
 print, 'MinTe=',MinTe, '   MaxTe=',MaxTe
 print, 'MinVsw=',MinVsw, '   MaxVsw=',MaxVsw
 print, 'Minkappa=',Minkappa, '   Maxkappa=',Maxkappa
;--------------------------------------------------------------------
      
;bin theta fluxes into 2D array: 
 theta=DBLARR(XH,2*YH-1)
 this_source=0L  
 temp[*,*]=0L                                
 FOR i=0,t_steps-1 DO BEGIN
  FOR j=0,Number_of_Sources[i]-1 DO BEGIN
   Xbin=LONG(CEIL((Sources[0,this_source]/Radial_Height[t_steps-1])*(XH-1)))
   IF(Sources[1,this_source] LT 0.0)THEN BEGIN
    Ybin=LONG(CEIL((Sources[1,this_source]/Radial_Height[t_steps-1])*(YH-1)+YH-1))
   ENDIF ELSE BEGIN 
    Ybin=LONG(FLOOR((Sources[1,this_source]/Radial_Height[t_steps-1])*(YH-1)+YH-1))
   ENDELSE 
   theta[Xbin,Ybin]=theta[Xbin,Ybin]+Sources[2,this_source] 
   temp[Xbin,Ybin]=temp[Xbin,Ybin] +1L
   this_source=this_source+1L                         
  ENDFOR ;Number_of_Sources
 ENDFOR ;t_steps  

 theta=theta/DOUBLE(temp+(temp EQ 0L))
;--------------------------------------------------------------------
   
;bin N_e fluxes into 2D array: 
 N_e=DBLARR(XH,2*YH-1)
 this_source=0L  
 temp[*,*]=0L                                
 FOR i=0,t_steps-1 DO BEGIN
  FOR j=0,Number_of_Sources[i]-1 DO BEGIN
   Xbin=LONG(CEIL((Sources[0,this_source]/Radial_Height[t_steps-1])*(XH-1)))
   IF(Sources[1,this_source] LT 0.0)THEN BEGIN
    Ybin=LONG(CEIL((Sources[1,this_source]/Radial_Height[t_steps-1])*(YH-1)+YH-1))
   ENDIF ELSE BEGIN 
    Ybin=LONG(FLOOR((Sources[1,this_source]/Radial_Height[t_steps-1])*(YH-1)+YH-1))
   ENDELSE 
   N_e[Xbin,Ybin]=N_e[Xbin,Ybin]+Sources[3,this_source] 
   temp[Xbin,Ybin]=temp[Xbin,Ybin] +1L
   this_source=this_source+1L                         
  ENDFOR ;Number_of_Sources
 ENDFOR ;t_steps  

 N_e=N_e/DOUBLE(temp+(temp EQ 0L))
;--------------------------------------------------------------------
 
;bin B1 fluxes into 2D array: 
 B1=DBLARR(XH,2*YH-1)
 this_source=0L   
 temp[*,*]=0L                               
 FOR i=0,t_steps-1 DO BEGIN
  FOR j=0,Number_of_Sources[i]-1 DO BEGIN
   Xbin=LONG(CEIL((Sources[0,this_source]/Radial_Height[t_steps-1])*(XH-1)))
   IF(Sources[1,this_source] LT 0.0)THEN BEGIN
    Ybin=LONG(CEIL((Sources[1,this_source]/Radial_Height[t_steps-1])*(YH-1)+YH-1))
   ENDIF ELSE BEGIN 
    Ybin=LONG(FLOOR((Sources[1,this_source]/Radial_Height[t_steps-1])*(YH-1)+YH-1))
   ENDELSE 
   B1[Xbin,Ybin]=B1[Xbin,Ybin]+Sources[4,this_source] 
   temp[Xbin,Ybin]=temp[Xbin,Ybin] +1L
   this_source=this_source+1L                         
  ENDFOR ;Number_of_Sources
 ENDFOR ;t_steps  

 B1=B1/DOUBLE(temp+(temp EQ 0L))
;--------------------------------------------------------------------
 
;bin Ti fluxes into 2D array: 
 Ti=DBLARR(XH,2*YH-1)
 this_source=0L      
 temp[*,*]=0L                             
 FOR i=0,t_steps-1 DO BEGIN
  FOR j=0,Number_of_Sources[i]-1 DO BEGIN
   Xbin=LONG(CEIL((Sources[0,this_source]/Radial_Height[t_steps-1])*(XH-1)))
   IF(Sources[1,this_source] LT 0.0)THEN BEGIN
    Ybin=LONG(CEIL((Sources[1,this_source]/Radial_Height[t_steps-1])*(YH-1)+YH-1))
   ENDIF ELSE BEGIN 
    Ybin=LONG(FLOOR((Sources[1,this_source]/Radial_Height[t_steps-1])*(YH-1)+YH-1))
   ENDELSE 
   Ti[Xbin,Ybin]=Ti[Xbin,Ybin]+Sources[5,this_source]
   temp[Xbin,Ybin]=temp[Xbin,Ybin] +1L 
   this_source=this_source+1L                         
  ENDFOR ;Number_of_Sources
 ENDFOR ;t_steps  

 Ti=Ti/DOUBLE(temp+(temp EQ 0L))
;--------------------------------------------------------------------
 
;bin Te fluxes into 2D array: 
 Te=DBLARR(XH,2*YH-1)
 this_source=0L    
 temp[*,*]=0L                              
 FOR i=0,t_steps-1 DO BEGIN
  FOR j=0,Number_of_Sources[i]-1 DO BEGIN
   Xbin=LONG(CEIL((Sources[0,this_source]/Radial_Height[t_steps-1])*(XH-1)))
   IF(Sources[1,this_source] LT 0.0)THEN BEGIN
    Ybin=LONG(CEIL((Sources[1,this_source]/Radial_Height[t_steps-1])*(YH-1)+YH-1))
   ENDIF ELSE BEGIN 
    Ybin=LONG(FLOOR((Sources[1,this_source]/Radial_Height[t_steps-1])*(YH-1)+YH-1))
   ENDELSE 
   Te[Xbin,Ybin]=Te[Xbin,Ybin]+Sources[6,this_source] 
   temp[Xbin,Ybin]=temp[Xbin,Ybin] +1L
   this_source=this_source+1L                         
  ENDFOR ;Number_of_Sources
 ENDFOR ;t_steps 

 Te=Te/DOUBLE(temp+(temp EQ 0L))
;--------------------------------------------------------------------
 
;bin Vsw fluxes into 2D array: 
 Vsw=DBLARR(XH,2*YH-1)
 this_source=0L  
 temp[*,*]=0L                       
 FOR i=0,t_steps-1 DO BEGIN
  FOR j=0,Number_of_Sources[i]-1 DO BEGIN
   Xbin=LONG(CEIL((Sources[0,this_source]/Radial_Height[t_steps-1])*(XH-1)))
   IF(Sources[1,this_source] LT 0.0)THEN BEGIN
    Ybin=LONG(CEIL((Sources[1,this_source]/Radial_Height[t_steps-1])*(YH-1)+YH-1))
   ENDIF ELSE BEGIN 
    Ybin=LONG(FLOOR((Sources[1,this_source]/Radial_Height[t_steps-1])*(YH-1)+YH-1))
   ENDELSE 
   Vsw[Xbin,Ybin]=Vsw[Xbin,Ybin]+Sources[7,this_source]
   temp[Xbin,Ybin]=temp[Xbin,Ybin] +1L 
   this_source=this_source+1L                         
  ENDFOR ;Number_of_Sources
 ENDFOR ;t_steps  

 Vsw=Vsw/DOUBLE(temp+(temp EQ 0L))
;--------------------------------------------------------------------
 
;bin kappa fluxes into 2D array: 
 kappa=DBLARR(XH,2*YH-1)
 this_source=0L   
 temp[*,*]=0L                               
 FOR i=0,t_steps-1 DO BEGIN
  FOR j=0,Number_of_Sources[i]-1 DO BEGIN
   Xbin=LONG(CEIL((Sources[0,this_source]/Radial_Height[t_steps-1])*(XH-1)))
   IF(Sources[1,this_source] LT 0.0)THEN BEGIN
    Ybin=LONG(CEIL((Sources[1,this_source]/Radial_Height[t_steps-1])*(YH-1)+YH-1))
   ENDIF ELSE BEGIN 
    Ybin=LONG(FLOOR((Sources[1,this_source]/Radial_Height[t_steps-1])*(YH-1)+YH-1))
   ENDELSE 
   kappa[Xbin,Ybin]=kappa[Xbin,Ybin]+Sources[8,this_source]
   temp[Xbin,Ybin]=temp[Xbin,Ybin] +1L
   this_source=this_source+1L                         
  ENDFOR ;Number_of_Sources
 ENDFOR ;t_steps 

 kappa=kappa/DOUBLE(temp+(temp EQ 0L))
;--------------------------------------------------------------------
  
;Obtain MAX & MIN nonzero [theta, N_e, B1, Ti, Te, Vsw, Kappa]: 
 Maxtheta=MAX(theta[WHERE(theta NE 0.0)],MIN=Mintheta) 
 MaxNe=MAX(N_e[WHERE(N_e NE 0.0)],MIN=MinNe) 
 MaxB1=MAX(B1[WHERE(B1 NE 0.0)],MIN=MinB1) 
 MaxTi=MAX(Ti[WHERE(Ti NE 0.0)],MIN=MinTi)   
 MaxTe=MAX(Te[WHERE(Te NE 0.0)],MIN=MinTe) 
 MaxVsw=MAX(Vsw[WHERE(Vsw NE 0.0)],MIN=MinVsw) 
 Maxkappa=MAX(kappa[WHERE(kappa NE 0.0)],MIN=Minkappa) 
 print, 'Mintheta=',Mintheta, '   Maxtheta=',Maxtheta
 print, 'MinNe=',MinNe, '   MaxNe=',MaxNe
 print, 'MinB1=',MinB1, '   MaxB1=',MaxB1
 print, 'MinTi=',MinTi, '   MaxTi=',MaxTi
 print, 'MinTe=',MinTe, '   MaxTe=',MaxTe
 print, 'MinVsw=',MinVsw, '   MaxVsw=',MaxVsw
 print, 'Minkappa=',Minkappa, '   Maxkappa=',Maxkappa

;set value=0.0 [theta, N_e, B1, Ti, Te, Vsw, Kappa] to minimum nonzero value: 
 theta=theta+(theta EQ 0.0)*Mintheta
 N_e=N_e+(N_e EQ 0.0)*MinNe
 B1=B1+(B1 EQ 0.0)*MinB1
 Ti=Ti+(Ti EQ 0.0)*MinTi
 Te=Te+(Te EQ 0.0)*MinTe
 Vsw=Vsw+(Vsw EQ 0.0)*MinVsw
 kappa=kappa+(kappa EQ 0.0)*Minkappa 
  
;take log of values:
 N_e=ALOG10(N_e)
 B1=ALOG10(B1)
 
;Obtain MAX & MIN in log form
 MaxNe=Max(N_e,MIN=MinNe)
 print, 'Log(MinNe)=',MinNe, '   Log(MaxNe)=',MaxNe
 MaxB1=Max(B1,MIN=MinB1)
 print, 'Log(MinB1)=',MinB1, '   Log(MaxB1)=',MaxB1

;get the range of values in  [theta, N_e, B1, Ti, Te, Vsw, Kappa]:
 theta_Range=Maxtheta-Mintheta
 Ne_Range=MaxNe-MinNe
 B1_Range=MaxB1-MinB1
 Ti_Range=MaxTi-MinTi
 Te_Range=MaxTe-MinTe
 Vsw_Range=MaxVsw-MinVsw
 Kappa_Range=MaxKappa-MinKappa
 print, 'Ne_Range=',Ne_Range
 print, 'B1_Range=',B1_Range
 print, 'Ti_Range=',Ti_Range
 print, 'Te_Range=',Te_Range
 print, 'Vsw_Range=',Vsw_Range
 print, 'Kappa_Range=',Kappa_Range
;--------------------------------------------------------------------

;array sizes and scaling factors
 scale=1.0e9
 XHs=350  &  YHs=350  
;--------------------------------------------------------------------

;Connect to the file containing the structure definitions
 openr, stru, outputDir+'parameters_'+fileident, /get_lun
  
 readf, stru, FORMAT='(49(/),10X,I2,4(/))', nos &  print, 'nos=',nos
 nos_string=STRTRIM(STRING(Fix(nos)),2)

 IF(nos GT 0)THEN BEGIN
   Structure_type=INTARR(nos) ;& print, 'Structure_type=',Structure_type
   par=FLTARR(nos,6) ;&  print, 'par=',par

   ;load the structures
   FOR i=0,nos-1 DO BEGIN
    readf, stru, FORMAT='(10X,I2,12X,I2)', Structure,temps  
    Structure_type[i]=temps
    CASE Structure_type[i] OF
     '1': begin
        readf, stru, FORMAT='(2(9X,F4.2),2(/))', temp1,temp2
        par[i,0]=temp1 & par[i,1]=temp2
     end
     '2': begin
        readf, stru, FORMAT='(/,8X,E12.3,/,8X,E12.3,3(/))', temp1,temp2
        par[i,0]=temp1 & par[i,1]=temp2
     end
     '3': begin
        readf, stru, FORMAT='(2X,E12.3,2(/),2X,E12.3,7X,E12.3,2(/))', $ 
               temp1,temp2,temp3
        par[i,0]=temp1 & par[i,1]=temp2 & par[i,2]=temp3
     end
     '4': begin
        readf, stru,                                                   $ 
               FORMAT='(/,6X,E12.3,11X,E12.3,/,6X,E12.3,11X,E12.3,2(/))', $ 
               temp1,temp2,temp3,temp4
        par[i,0]=temp1 & par[i,1]=temp2 & par[i,2]=temp3 & par[i,3]=temp4  
     end
     '6': begin
        readf, stru, FORMAT=                                                  $ 
          '(/,3X,E12.3,/,3X,E12.3,8X,E12.3,/,3X,E12.3,/,3X,E12.3,8X,E12.3,2(/))', $ 
               temp1,temp2,temp3,temp4,temp5,temp6 
        par[i,0]=temp1 & par[i,1]=temp2 & par[i,2]=temp3 
        par[i,3]=temp4 & par[i,4]=temp5 & par[i,5]=temp6 
     end
    ELSE: begin
            print, 'seems that the structure type is not recognised'
            stop
          end
    ENDCASE
   ENDFOR ;load the structures
   free_lun,stru
   print, 'par=',par 
  ;--------------------------------------------------------------------

  ;Set up arrays 
   X=DBLARR(XHs,2*YHs-1) 
   FOR j=LONG(0),XHs-1 DO BEGIN 
     X[j,*]=j+1
   ENDFOR ;
   X=X*RangeGm/DOUBLE(XHs)
   Y=DBLARR(XHs,2*YHs-1) 
   FOR k=LONG(0),2*YHs-2 DO BEGIN 
     Y[*,k]=k-YHs+1
   ENDFOR ;
   Y=Y*RangeGm/DOUBLE(YHs)
   RH=SQRT(X*X + Y*Y)
   R_spiral=137*(0.5*pi-ATAN(Y/X))*Sol_rad/scale
   Structures=INTARR(XHs,2*YHs-1,nos) 
  ;--------------------------------------------------------------------

  ;Put Structures into arrays
   FOR i=0,nos-1 DO BEGIN 
     CASE Structure_type[i] OF
      '1': begin  ;CIR
         R_spiral=137*(par[i,0]*pi-ATAN(Y/X))*Sol_rad/scale
         Index=WHERE((((R_spiral-RH) GT 0) AND ((R_spiral*(1-par[i,1])-RH) LT 0)) $ 
                OR (((R_spiral-RH) LT 0) AND ((R_spiral*(1+par[i,1])-RH) GT 0))   $ 
                OR (((R_spiral*(1-par[i,1])-RH) GT 0)                             $ 
                AND ((R_spiral*(1-2*par[i,1])-RH) LT 0)))
         index_size=SIZE(index)
         IF(index_size[0] eq 1)THEN BEGIN   
          Yvals=LONG(Index/DOUBLE(XHs))   
          Xvals=Index-Yvals*XHs
          FOR jj=LONG(0),index_size[1]-1 DO BEGIN     
           Structures[Xvals[jj],Yvals[jj],i]=i+1
          ENDFOR
         ENDIF ELSE print, 'Failed to resolve structure ',i,' which is a type 1'
      end
      '2': begin  ;annulus between heliocentric radius aaa & bbb
         Index=WHERE(((RH*RH) GT (par[i,0]*par[i,0]/(scale*scale))) $ 
                     AND ((RH*RH) LT (par[i,1]*par[i,1]/(scale*scale))))
         index_size=SIZE(index)
         IF(index_size[0] eq 1)THEN BEGIN   
          Yvals=LONG(Index/DOUBLE(XHs))   
          Xvals=Index-Yvals*XHs
          FOR jj=LONG(0),index_size[1]-1 DO BEGIN     
           Structures[Xvals[jj],Yvals[jj],i]=i+1
          ENDFOR
         ENDIF ELSE print, 'Failed to resolve structure ',i,' which is a type 2'
      end
      '4': begin   ;square/rectangular region bound by aaa, bbb, ccc, ddd
         Index=WHERE((X GT par[i,0]/scale) AND (X LT par[i,1]/scale) $ 
                    AND (Y LT par[i,2]/scale) AND (Y GT par[i,3]/scale)) 
         index_size=SIZE(index) 
         IF(index_size[0] eq 1)THEN BEGIN     
          Yvals=LONG(Index/DOUBLE(XHs))   
          Xvals=Index-Yvals*XHs
          FOR jj=LONG(0),index_size[1]-1 DO BEGIN   
           Structures[Xvals[jj],Yvals[jj],i]=i+1
          ENDFOR
         ENDIF ELSE print, 'Failed to resolve structure ',i,' which is a type 4'
      end 
      '3': begin ;circular region of radius aaa, & center (X,Y)=(bbb,ccc)
         Index=WHERE(((X-par[i,1]/scale)^2+(Y-par[i,2]/scale)^2) $ 
                      LT ((par[i,0]*par[i,0])/(scale*scale)))
         index_size=SIZE(index)
         IF(index_size[0] eq 1)THEN BEGIN     
          Yvals=LONG(Index/DOUBLE(XHs))   
          Xvals=Index-Yvals*XHs
          FOR jj=LONG(0),index_size[1]-1 DO BEGIN   
           Structures[Xvals[jj],Yvals[jj],i]=i+1
          ENDFOR
         ENDIF ELSE print, 'Failed to resolve structure ',i,' which is a type 3'
      end 
      '6': begin ;coronal loop
         Index=WHERE(((X-par[i,1]/scale)^2+(Y-par[i,2]/scale)^2)     $ 
                     GT ((par[i,0]*par[i,0])/(scale*scale))          $
                     AND ((X-par[i,4]/scale)^2+(Y-par[i,5]/scale)^2) $
                     LT ((par[i,3]*par[i,3])/(scale*scale)))
         index_size=SIZE(index)
         IF(index_size[0] eq 1)THEN BEGIN      
           Yvals=LONG(Index/DOUBLE(XHs))   
           Xvals=Index-Yvals*XHs
           FOR jj=LONG(0),index_size[1]-1 DO BEGIN   
             Structures[Xvals[jj],Yvals[jj],i]=i+1
           ENDFOR
         ENDIF ELSE print, 'Failed to resolve structure ',i,' which is a type 6'
          end 
     ELSE: begin
             print, 'something has gone wrong with the structure type...'
             STOP
           end
     ENDCASE
   ENDFOR ;Put Structures into arrays
 ENDIF
;-------------------------------------------------------------------- 

;read in all three observer locations 
 openr, lun, outputDir+'parameters_'+fileident, /get_lun
 readf, lun,  FORMAT= $ 
  '(32(/),2X,E11.2,/,2X,E12.2,3(/),2X,E11.2,/,2X,E12.2,3(/),2X,E11.2,/,2X,E12.2)', $ 
   Obs1X, Obs1Y, Obs2X, Obs2Y, Obs3X, Obs3Y
 free_lun,lun

 CASE Observer OF
   '1': begin
          ObserverX=Obs1X/xscale
          ObserverY=Obs1Y/yscale 
        end
   '2': begin
          ObserverX=Obs2X/xscale
          ObserverY=Obs2Y/yscale 
        end
   '3': begin
          ObserverX=Obs3X/xscale
          ObserverY=Obs3Y/yscale 
        end
 ELSE: begin
         print, 'That Observer does not exist, sorry.'
         stop
       end
 ENDCASE
 print, ObserverX, ObserverY
;--------------------------------------------------------------------

;bin Fund fluxes into 2D array: 
 Fund=DBLARR(XH,2*YH-1)
 this_source=0L                         
 FOR i=0,t_steps-1 DO BEGIN
 FOR j=0,Number_of_Sources[i]-1 DO BEGIN
  Xbin=LONG(CEIL((SourcesFH[0,this_source]/Radial_Height[t_steps-1])*(XH-1)))
  IF(Sources[1,this_source] LT 0.0)THEN BEGIN
   Ybin=LONG(CEIL((SourcesFH[1,this_source]/Radial_Height[t_steps-1])*(YH-1)+YH-1))
  ENDIF ELSE BEGIN 
   Ybin=LONG(FLOOR((SourcesFH[1,this_source]/Radial_Height[t_steps-1])*(YH-1)+YH-1))
  ENDELSE 
  Fund[Xbin,Ybin]=Fund[Xbin,Ybin]+SourcesFH[2,this_source] 
  this_source=this_source+1L                         
 ENDFOR ;Number_of_Sources
 ENDFOR ;t_steps 
;--------------------------------------------------------------------
  
;bin harmonic fluxes into 2D array: 
 Harm=DBLARR(XH,2*YH-1)
 this_source=0L                         
 FOR i=0,t_steps-1 DO BEGIN
 FOR j=0,Number_of_Sources[i]-1 DO BEGIN
  Xbin=LONG(CEIL((SourcesFH[0,this_source]/Radial_Height[t_steps-1])*(XH-1)))
  IF(Sources[1,this_source] LT 0.0)THEN BEGIN
   Ybin=LONG(CEIL((SourcesFH[1,this_source]/Radial_Height[t_steps-1])*(YH-1)+YH-1))
  ENDIF ELSE BEGIN 
   Ybin=LONG(FLOOR((SourcesFH[1,this_source]/Radial_Height[t_steps-1])*(YH-1)+YH-1))
  ENDELSE 
  Harm[Xbin,Ybin]=Harm[Xbin,Ybin]+SourcesFH[3,this_source] 
  this_source=this_source+1L                         
 ENDFOR ;Number_of_Sources
 ENDFOR ;t_steps 
;--------------------------------------------------------------------
   
;Obtain MAX & MIN nonzero Fund & Harm Fluxes: 
 MaxF=MAX(Fund[WHERE(Fund NE 0.0)],MIN=MinF) 
 MaxH=MAX(Harm[WHERE(Harm NE 0.0)],MIN=MinH) 
 print, 'MinF=',MinF, '   MaxF=',MaxF
 print, 'MinH=',MinH, '   MaxH=',MaxH
;--------------------------------------------------------------------

;set value=0.0 fund & harm fluxes to minimum nonzero value: 
 Fund=Fund+(Fund EQ 0.0)*MinF
 Harm=Harm+(Harm EQ 0.0)*MinH 
;-------------------------------------------------------------------- 
  
;take log of fluxes so that values are ~ -21.0 instead of 1x10^(-21.0):
 Fund=ALOG10(Fund)
 Harm=ALOG10(Harm)
;--------------------------------------------------------------------
 
;Obtain MAX & MIN in log form
 MaxF=Max(Fund,MIN=MinF)
 MaxH=Max(Harm,MIN=MinH) 
 print, 'MinF=',MinF, '   MaxF=',MaxF
 print, 'MinH=',MinH, '   MaxH=',MaxH
;--------------------------------------------------------------------

;get min & max of flux range that contains both Fund & Harm
 MinHF=MAX([MinF, MinH]) 
 MaxHF=MAX([MaxF, MaxH]) 
 print, 'MinHF=',MinHF
 print, 'MaxHF=',MaxHF
;--------------------------------------------------------------------
 
;Set lower limit of flux for plots:
 Fund=Fund > (MinHF)
 Harm=Harm > (MinHF)
;--------------------------------------------------------------------

;get the range of flux values in fund & harm:
 HF_Range=MaxHF-MinHF
 print, 'HF_Range=',HF_Range
;--------------------------------------------------------------------
;____________________________________________________________________
  

;--------------------------------------------------------------------
;******************* Initialize output device ***********************
;====================================================================
               
 OPFilename='Helio'+observer
 NHeaderRows=0
 x_cm=27.0   ; < 21.0 or 29.7
 y_cm=18.7   ; < 29.7 or 21.0
 NCLev=7 
 NGLev=8*NCLev

;set plot device:                 
 SelectDevice, OPFilename, NHeaderRows, x_cm, y_cm, OPDevice, NCLev, NGLev
 IF(OPDevice ne 'q')THEN  begin
   !P.Multi=[0,5,2]
    
  ;OVER-RIDE  object margins
   !x.omargin=[0,0]     ; [left,right]
   !y.omargin=[1,0]     ; [bottom,top]

  ;OVER-RIDE  plot margins
   !x.margin=[8,4]     ; [left,right]
   !y.margin=[7,3]     ; [bottom,top]
        
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
        
   restore, 'colours/DS_WG.sav'
   tvlct, r, g, b            
  ;--------------------------------------------------------------------
  ;____________________________________________________________________


  ;--------------------------------------------------------------------
  ;******************************** Plot ******************************
  ;==================================================================== 

  ;set levels for contour plot
   stepg=theta_Range/DOUBLE(NGLev);
   grey=indgen(NGLev)*stepg+Mintheta;

  ;theta 
   contour, theta, (FIndGen(XH)*Radial_Height[t_steps-1])/DOUBLE(XH-1)/xscale,   $
                  ((FIndGen(2*YH-1)*Radial_Height[t_steps-1])/DOUBLE(YH-1)       $ 
                   -Radial_Height[t_steps-1])/yscale,                            $
                   Title=TITLE_FONT+' theta !3 (degrees)',                       $ 
                   XTitle='X!DH!N  (Gm)', YTitle='Y!DH!N  (Gm)', charsize=1.7,   $
                   XRange=[0,RangeGm], YRange=[-RangeGm,RangeGm],                $
                   levels=grey, /fill, c_colors=1+indgen(NGLev)

  ;Global shocks at periodic time steps
   FOR i=0,t_steps-1, k DO BEGIN
   IF((Radial_Height[i] GT sol_rad) AND (Radial_Height[i] LT RangeGm*1.0e9))THEN BEGIN
     x_shk=FINDGEN(arr)*(Radial_Height[i])/DOUBLE(arr-1)
     Y_shk=SQRT(x_shk/CMEb[i]) & Y_shk= ROTATE(Y_shk,2)
     Y2_shk= -Y_shk
     j=LONG(sol_rad/(Radial_Height[i]/DOUBLE(arr-1)))
     Plots, x_shk[j:(arr-1)]/xscale, Y_shk[j:(arr-1)]/yscale, THICK=3 
     Plots, x_shk[j:(arr-1)]/xscale, Y2_shk[j:(arr-1)]/yscale, THICK=3
   ENDIF
   ENDFOR

  ;Sun
   x=FINDGEN(arr)*(sol_rad/xscale)/DOUBLE(arr-1)
   OPLOT, x, SQRT((sol_rad/xscale)^2-x^2),  color=254, THICK=4
   OPLOT, x, -SQRT((sol_rad/xscale)^2-x^2), color=254, THICK=4

  ;Concentric circles at integer numbers of solar radii
   FOR i=2,5 DO BEGIN
     x=FINDGEN(arr)*(i*sol_rad/xscale)/DOUBLE(arr-1)
     OPLOT, x, SQRT((i*sol_rad/xscale)^2-x^2),  color=254, linestyle=1 
     OPLOT, x, -SQRT((i*sol_rad/xscale)^2-x^2), color=254, linestyle=1 
   ENDFOR
   FOR i=10,25, 5 DO BEGIN
     x=FINDGEN(arr)*(i*sol_rad/xscale)/DOUBLE(arr-1)
     OPLOT, x, SQRT((i*sol_rad/xscale)^2-x^2),  color=254, linestyle=1 
     OPLOT, x, -SQRT((i*sol_rad/xscale)^2-x^2), color=254, linestyle=1 
   ENDFOR
   FOR i=50,200, 50 DO BEGIN
     x=FINDGEN(arr)*(i*sol_rad/xscale)/DOUBLE(arr-1)
     OPLOT, x, SQRT((i*sol_rad/xscale)^2-x^2),  color=254, linestyle=1 
     OPLOT, x, -SQRT((i*sol_rad/xscale)^2-x^2), color=254, linestyle=1 
   ENDFOR

   ;______________
   ;
   ; Add colourbar 
   ;______________ 

   COLORBAR, CHARSIZE = 1.5,  DIVISIONS = 5, FORMAT = '(I3)', $ 
   POSITION = [0.03, 0.519, 0.18, 0.539],  $ 
   NCOLORS = NGLev-1, RANGE = [180L*mintheta/Pi , 180L*maxtheta/Pi ], BOTTOM = 2
  ;--------------------------------------------------------------------

  ;set levels for contour plot
   stepg=Ne_Range/DOUBLE(NGLev);
   grey=indgen(NGLev)*stepg+MinNe;

  ;N_e Emission intensities
   contour, N_e, (FIndGen(XH)*Radial_Height[t_steps-1])/DOUBLE(XH-1)/xscale,  $ 
                 ((FIndGen(2*YH-1)*Radial_Height[t_steps-1])/DOUBLE(YH-1)     $ 
                 -Radial_Height[t_steps-1])/yscale,                           $
                 Title=TITLE_FONT+' N!De!N !3  Log(m!U-3!N)',                 $
                 XTitle='X!DH!N  (Gm)', YTitle='Y!DH!N  (Gm)', charsize=1.7,  $ 
                 XRange=[0,RangeGm], YRange=[-RangeGm,RangeGm],               $
                 levels=grey, /fill, c_colors=1+indgen(NGLev)

  ;Global shocks at periodic time steps
   FOR i=0,t_steps-1, k DO BEGIN
   IF((Radial_Height[i] GT sol_rad) AND (Radial_Height[i] LT RangeGm*1.0e9))THEN BEGIN
    x_shk=FINDGEN(arr)*Radial_Height[i]/DOUBLE(arr-1)
    Y_shk=SQRT(x_shk/CMEb[i]) & Y_shk= ROTATE(Y_shk,2)
    Y2_shk= -Y_shk
    j=LONG(sol_rad/(Radial_Height[i]/DOUBLE(arr-1)))   
    Plots, x_shk[j:(arr-1)]/xscale, Y_shk[j:(arr-1)]/yscale, THICK=3 
    Plots, x_shk[j:(arr-1)]/xscale, Y2_shk[j:(arr-1)]/yscale, THICK=3
   ENDIF
   ENDFOR

  ;Sun
   x=FINDGEN(arr)*(sol_rad/xscale)/DOUBLE(arr-1)
   OPLOT, x, SQRT((sol_rad/xscale)^2-x^2),  color=254, THICK=4
   OPLOT, x, -SQRT((sol_rad/xscale)^2-x^2), color=254, THICK=4

  ;Concentric circles at integer numbers of solar radii
   FOR i=2,5 DO BEGIN
     x=FINDGEN(arr)*(i*sol_rad/xscale)/DOUBLE(arr-1)
     OPLOT, x, SQRT((i*sol_rad/xscale)^2-x^2),  color=254, linestyle=1 
     OPLOT, x, -SQRT((i*sol_rad/xscale)^2-x^2), color=254, linestyle=1 
   ENDFOR
   FOR i=10,25, 5 DO BEGIN
     x=FINDGEN(arr)*(i*sol_rad/xscale)/DOUBLE(arr-1)
     OPLOT, x, SQRT((i*sol_rad/xscale)^2-x^2),  color=254, linestyle=1 
     OPLOT, x, -SQRT((i*sol_rad/xscale)^2-x^2), color=254, linestyle=1 
   ENDFOR
   FOR i=50,200, 50 DO BEGIN
     x=FINDGEN(arr)*(i*sol_rad/xscale)/DOUBLE(arr-1)
     OPLOT, x, SQRT((i*sol_rad/xscale)^2-x^2),  color=254, linestyle=1 
     OPLOT, x, -SQRT((i*sol_rad/xscale)^2-x^2), color=254, linestyle=1 
   ENDFOR

   ;______________
   ;
   ; Add colourbar 
   ;______________ 

   COLORBAR,  CHARSIZE = 1.5, DIVISIONS = 5, FORMAT = '(F5.2)', $ 
   POSITION = [0.23, 0.519, 0.38, 0.539],  $ 
   NCOLORS = NGLev-1, RANGE = [minNe, maxNe], BOTTOM = 2
  ;--------------------------------------------------------------------

  ;set levels for contour plot
   stepg=B1_Range/DOUBLE(NGLev);
   grey=indgen(NGLev)*stepg+MinB1;


  ;B1 Emission intensities
   contour, B1, (FIndGen(XH)*Radial_Height[t_steps-1])/DOUBLE(XH-1)/xscale,  $ 
                ((FIndGen(2*YH-1)*Radial_Height[t_steps-1])/DOUBLE(YH-1)     $ 
                -Radial_Height[t_steps-1])/yscale,                           $
                Title=TITLE_FONT+' B !3  Log(T)',                            $
                XTitle='X!DH!N  (Gm)', YTitle='Y!DH!N  (Gm)', charsize=1.7,  $ 
                XRange=[0,RangeGm], YRange=[-RangeGm,RangeGm],               $
                levels=grey, /fill, c_colors=1+indgen(NGLev)

  ;Global shocks at periodic time steps
   FOR i=0,t_steps-1, k DO BEGIN
   IF((Radial_Height[i] GT sol_rad) AND (Radial_Height[i] LT RangeGm*1.0e9))THEN BEGIN
    x_shk=FINDGEN(arr)*Radial_Height[i]/DOUBLE(arr-1)
    Y_shk=SQRT(x_shk/CMEb[i]) & Y_shk= ROTATE(Y_shk,2)
    Y2_shk= -Y_shk
    j=LONG(sol_rad/(Radial_Height[i]/DOUBLE(arr-1)))   
    Plots, x_shk[j:(arr-1)]/xscale, Y_shk[j:(arr-1)]/yscale, THICK=3 
    Plots, x_shk[j:(arr-1)]/xscale, Y2_shk[j:(arr-1)]/yscale, THICK=3
   ENDIF
   ENDFOR

  ;Sun
   x=FINDGEN(arr)*(sol_rad/xscale)/DOUBLE(arr-1)
   OPLOT, x, SQRT((sol_rad/xscale)^2-x^2),  color=254, THICK=4
   OPLOT, x, -SQRT((sol_rad/xscale)^2-x^2), color=254, THICK=4 

  ;Concentric circles at integer numbers of solar radii
   FOR i=2,5 DO BEGIN
     x=FINDGEN(arr)*(i*sol_rad/xscale)/DOUBLE(arr-1)
     OPLOT, x, SQRT((i*sol_rad/xscale)^2-x^2),  color=254, linestyle=1 
     OPLOT, x, -SQRT((i*sol_rad/xscale)^2-x^2), color=254, linestyle=1 
   ENDFOR
   FOR i=10,25, 5 DO BEGIN
     x=FINDGEN(arr)*(i*sol_rad/xscale)/DOUBLE(arr-1)
     OPLOT, x, SQRT((i*sol_rad/xscale)^2-x^2),  color=254, linestyle=1 
     OPLOT, x, -SQRT((i*sol_rad/xscale)^2-x^2), color=254, linestyle=1 
   ENDFOR
   FOR i=50,200, 50 DO BEGIN
     x=FINDGEN(arr)*(i*sol_rad/xscale)/DOUBLE(arr-1)
     OPLOT, x, SQRT((i*sol_rad/xscale)^2-x^2),  color=254, linestyle=1 
     OPLOT, x, -SQRT((i*sol_rad/xscale)^2-x^2), color=254, linestyle=1 
   ENDFOR

   ;______________
   ;
   ; Add colourbar 
   ;______________ 

   COLORBAR,  CHARSIZE = 1.5, DIVISIONS = 5, FORMAT = '(F4.1)', $ 
   POSITION = [0.43, 0.519, 0.58, 0.539],  $ 
   NCOLORS = NGLev-1, RANGE = [minB1, maxB1], BOTTOM = 2
  ;--------------------------------------------------------------------

  ;set levels for contour plot
   stepg=Ti_Range/DOUBLE(NGLev);
   grey=indgen(NGLev)*stepg+MinTi;

  ;Ti Emission intensities
   contour, Ti, (FIndGen(XH)*Radial_Height[t_steps-1])/DOUBLE(XH-1)/xscale,  $ 
                ((FIndGen(2*YH-1)*Radial_Height[t_steps-1])/DOUBLE(YH-1)     $ 
                -Radial_Height[t_steps-1])/yscale,                           $
                Title=TITLE_FONT+' T!Di!N !3 (10!U4!N K)',                   $
                XTitle='X!DH!N  (Gm)', YTitle='Y!DH!N  (Gm)', charsize=1.7,  $ 
                XRange=[0,RangeGm], YRange=[-RangeGm,RangeGm],               $
                levels=grey, /fill, c_colors=1+indgen(NGLev)

  ;Global shocks at periodic time steps
   FOR i=0,t_steps-1, k DO BEGIN
   IF((Radial_Height[i] GT sol_rad) AND (Radial_Height[i] LT RangeGm*1.0e9))THEN BEGIN
    x_shk=FINDGEN(arr)*Radial_Height[i]/DOUBLE(arr-1)
    Y_shk=SQRT(x_shk/CMEb[i]) & Y_shk= ROTATE(Y_shk,2)
    Y2_shk= -Y_shk
    j=LONG(sol_rad/(Radial_Height[i]/DOUBLE(arr-1)))   
    Plots, x_shk[j:(arr-1)]/xscale, Y_shk[j:(arr-1)]/yscale, THICK=3 
    Plots, x_shk[j:(arr-1)]/xscale, Y2_shk[j:(arr-1)]/yscale, THICK=3
   ENDIF
   ENDFOR

  ;Sun
   x=FINDGEN(arr)*(sol_rad/xscale)/DOUBLE(arr-1)
   OPLOT, x, SQRT((sol_rad/xscale)^2-x^2),  color=254, THICK=4
   OPLOT, x, -SQRT((sol_rad/xscale)^2-x^2), color=254, THICK=4 

  ;Concentric circles at integer numbers of solar radii
   FOR i=2,5 DO BEGIN
     x=FINDGEN(arr)*(i*sol_rad/xscale)/DOUBLE(arr-1)
     OPLOT, x, SQRT((i*sol_rad/xscale)^2-x^2),  color=254, linestyle=1 
     OPLOT, x, -SQRT((i*sol_rad/xscale)^2-x^2), color=254, linestyle=1 
   ENDFOR
   FOR i=10,25, 5 DO BEGIN
     x=FINDGEN(arr)*(i*sol_rad/xscale)/DOUBLE(arr-1)
     OPLOT, x, SQRT((i*sol_rad/xscale)^2-x^2),  color=254, linestyle=1 
     OPLOT, x, -SQRT((i*sol_rad/xscale)^2-x^2), color=254, linestyle=1 
   ENDFOR
   FOR i=50,200, 50 DO BEGIN
     x=FINDGEN(arr)*(i*sol_rad/xscale)/DOUBLE(arr-1)
     OPLOT, x, SQRT((i*sol_rad/xscale)^2-x^2),  color=254, linestyle=1 
     OPLOT, x, -SQRT((i*sol_rad/xscale)^2-x^2), color=254, linestyle=1 
   ENDFOR

   ;______________
   ;
   ; Add colourbar 
   ;______________ 

   COLORBAR,  CHARSIZE = 1.5, DIVISIONS = 5, FORMAT = '(I5)', $ 
   POSITION = [0.63, 0.519, 0.78, 0.539],  $ 
   NCOLORS = NGLev-1, RANGE = [minTi/10000.0D, maxTi/10000.0D], BOTTOM = 2
  ;--------------------------------------------------------------------

  ;set levels for contour plot
   stepg=Te_Range/DOUBLE(NGLev);
   grey=indgen(NGLev)*stepg+MinTe;

  ;Te Emission intensities
   contour, Te, (FIndGen(XH)*Radial_Height[t_steps-1])/DOUBLE(XH-1)/xscale,  $ 
                ((FIndGen(2*YH-1)*Radial_Height[t_steps-1])/DOUBLE(YH-1)     $ 
                -Radial_Height[t_steps-1])/yscale,                           $
                Title=TITLE_FONT+' T!De!N !3 (10!U4!N K)',                   $
                XTitle='X!DH!N  (Gm)', YTitle='Y!DH!N  (Gm)', charsize=1.7,  $ 
                XRange=[0,RangeGm], YRange=[-RangeGm,RangeGm],               $
                levels=grey, /fill, c_colors=1+indgen(NGLev)

  ;Global shocks at periodic time steps
   FOR i=0,t_steps-1, k DO BEGIN
   IF((Radial_Height[i] GT sol_rad) AND (Radial_Height[i] LT RangeGm*1.0e9))THEN BEGIN
    x_shk=FINDGEN(arr)*Radial_Height[i]/DOUBLE(arr-1)
    Y_shk=SQRT(x_shk/CMEb[i]) & Y_shk= ROTATE(Y_shk,2)
    Y2_shk= -Y_shk
    j=LONG(sol_rad/(Radial_Height[i]/DOUBLE(arr-1)))   
    Plots, x_shk[j:(arr-1)]/xscale, Y_shk[j:(arr-1)]/yscale, THICK=3 
    Plots, x_shk[j:(arr-1)]/xscale, Y2_shk[j:(arr-1)]/yscale, THICK=3
   ENDIF
   ENDFOR

  ;Sun
   x=FINDGEN(arr)*(sol_rad/xscale)/DOUBLE(arr-1)
   OPLOT, x, SQRT((sol_rad/xscale)^2-x^2),  color=254, THICK=4
   OPLOT, x, -SQRT((sol_rad/xscale)^2-x^2), color=254, THICK=4 

  ;Concentric circles at integer numbers of solar radii
   FOR i=2,5 DO BEGIN
     x=FINDGEN(arr)*(i*sol_rad/xscale)/DOUBLE(arr-1)
     OPLOT, x, SQRT((i*sol_rad/xscale)^2-x^2),  color=254, linestyle=1 
     OPLOT, x, -SQRT((i*sol_rad/xscale)^2-x^2), color=254, linestyle=1 
   ENDFOR
   FOR i=10,25, 5 DO BEGIN
     x=FINDGEN(arr)*(i*sol_rad/xscale)/DOUBLE(arr-1)
     OPLOT, x, SQRT((i*sol_rad/xscale)^2-x^2),  color=254, linestyle=1 
     OPLOT, x, -SQRT((i*sol_rad/xscale)^2-x^2), color=254, linestyle=1 
   ENDFOR
   FOR i=50,200, 50 DO BEGIN
     x=FINDGEN(arr)*(i*sol_rad/xscale)/DOUBLE(arr-1)
     OPLOT, x, SQRT((i*sol_rad/xscale)^2-x^2),  color=254, linestyle=1 
     OPLOT, x, -SQRT((i*sol_rad/xscale)^2-x^2), color=254, linestyle=1 
   ENDFOR

   ;______________
   ;
   ; Add colourbar 
   ;______________ 

   COLORBAR,  CHARSIZE = 1.5, DIVISIONS = 5, FORMAT = '(I4)', $ 
   POSITION = [0.83, 0.519, 0.98, 0.539],  $  
   NCOLORS = NGLev-1, RANGE = [minTe/10000.0D, maxTe/10000.0D], BOTTOM = 2
  ;--------------------------------------------------------------------

  ;set levels for contour plot
   stepg=Vsw_Range/DOUBLE(NGLev);
   grey=indgen(NGLev)*stepg+MinVsw;

  ;Vsw Emission intensities
   contour, Vsw, (FIndGen(XH)*Radial_Height[t_steps-1])/DOUBLE(XH-1)/xscale,   $ 
                 ((FIndGen(2*YH-1)*Radial_Height[t_steps-1])/DOUBLE(YH-1)      $ 
                 -Radial_Height[t_steps-1])/yscale,                            $
                 Title=TITLE_FONT+' Vsw !3 (km s!U-1!N)',                      $
                 XTitle='X!DH!N  (Gm)', YTitle='Y!DH!N  (Gm)', charsize=1.7,   $ 
                 XRange=[0,RangeGm], YRange=[-RangeGm,RangeGm],                $
                 levels=grey, /fill, c_colors=1+indgen(NGLev)

  ;Global shocks at periodic time steps
   FOR i=0,t_steps-1, k DO BEGIN
   IF((Radial_Height[i] GT sol_rad) AND (Radial_Height[i] LT RangeGm*1.0e9))THEN BEGIN
    x_shk=FINDGEN(arr)*Radial_Height[i]/DOUBLE(arr-1)
    Y_shk=SQRT(x_shk/CMEb[i]) & Y_shk= ROTATE(Y_shk,2)
    Y2_shk= -Y_shk
    j=LONG(sol_rad/(Radial_Height[i]/DOUBLE(arr-1)))   
    Plots, x_shk[j:(arr-1)]/xscale, Y_shk[j:(arr-1)]/yscale, THICK=3 
    Plots, x_shk[j:(arr-1)]/xscale, Y2_shk[j:(arr-1)]/yscale, THICK=3
   ENDIF
   ENDFOR

  ;Sun
   x=FINDGEN(arr)*(sol_rad/xscale)/DOUBLE(arr-1)
   OPLOT, x, SQRT((sol_rad/xscale)^2-x^2),  color=254, THICK=4
   OPLOT, x, -SQRT((sol_rad/xscale)^2-x^2), color=254, THICK=4 

  ;Concentric circles at integer numbers of solar radii
   FOR i=2,5 DO BEGIN
     x=FINDGEN(arr)*(i*sol_rad/xscale)/DOUBLE(arr-1)
     OPLOT, x, SQRT((i*sol_rad/xscale)^2-x^2),  color=254, linestyle=1 
     OPLOT, x, -SQRT((i*sol_rad/xscale)^2-x^2), color=254, linestyle=1 
   ENDFOR
   FOR i=10,25, 5 DO BEGIN
     x=FINDGEN(arr)*(i*sol_rad/xscale)/DOUBLE(arr-1)
     OPLOT, x, SQRT((i*sol_rad/xscale)^2-x^2),  color=254, linestyle=1 
     OPLOT, x, -SQRT((i*sol_rad/xscale)^2-x^2), color=254, linestyle=1 
   ENDFOR
   FOR i=50,200, 50 DO BEGIN
     x=FINDGEN(arr)*(i*sol_rad/xscale)/DOUBLE(arr-1)
     OPLOT, x, SQRT((i*sol_rad/xscale)^2-x^2),  color=254, linestyle=1 
     OPLOT, x, -SQRT((i*sol_rad/xscale)^2-x^2), color=254, linestyle=1 
   ENDFOR

   ;______________
   ;
   ; Add colourbar 
   ;______________ 

   COLORBAR,  CHARSIZE = 1.5, DIVISIONS = 5, FORMAT = '(I4)', $ 
   POSITION = [0.03, 0.023, 0.18, 0.043],  $ 
   NCOLORS = NGLev-1, RANGE = [minVsw/1000.0D, maxVsw/1000.0D], BOTTOM = 2
  ;--------------------------------------------------------------------

  ;set levels for contour plot
   stepg=kappa_Range/DOUBLE(NGLev);
   grey=indgen(NGLev)*stepg+Minkappa;

  ;kappa Emission intensities
   contour, kappa, (FIndGen(XH)*Radial_Height[t_steps-1])/DOUBLE(XH-1)/xscale,  $ 
                   ((FIndGen(2*YH-1)*Radial_Height[t_steps-1])/DOUBLE(YH-1)     $ 
                   -Radial_Height[t_steps-1])/yscale,                           $
                   Title=TITLE_FONT+' kappa !3',                                $
                   XTitle='X!DH!N  (Gm)', YTitle='Y!DH!N  (Gm)', charsize=1.7,  $ 
                   XRange=[0,RangeGm], YRange=[-RangeGm,RangeGm],               $
                   levels=grey, /fill, c_colors=1+indgen(NGLev)

  ;Global shocks at periodic time steps
   FOR i=0,t_steps-1, k DO BEGIN
   IF((Radial_Height[i] GT sol_rad) AND (Radial_Height[i] LT RangeGm*1.0e9))THEN BEGIN
    x_shk=FINDGEN(arr)*Radial_Height[i]/DOUBLE(arr-1)
    Y_shk=SQRT(x_shk/CMEb[i]) & Y_shk= ROTATE(Y_shk,2)
    Y2_shk= -Y_shk
    j=LONG(sol_rad/(Radial_Height[i]/DOUBLE(arr-1)))   
    Plots, x_shk[j:(arr-1)]/xscale, Y_shk[j:(arr-1)]/yscale, THICK=3 
    Plots, x_shk[j:(arr-1)]/xscale, Y2_shk[j:(arr-1)]/yscale, THICK=3
   ENDIF
   ENDFOR

  ;Sun
   x=FINDGEN(arr)*(sol_rad/xscale)/DOUBLE(arr-1)
   OPLOT, x, SQRT((sol_rad/xscale)^2-x^2),  color=254, THICK=4
   OPLOT, x, -SQRT((sol_rad/xscale)^2-x^2), color=254, THICK=4

  ;Concentric circles at integer numbers of solar radii
   FOR i=2,5 DO BEGIN
     x=FINDGEN(arr)*(i*sol_rad/xscale)/DOUBLE(arr-1)
     OPLOT, x, SQRT((i*sol_rad/xscale)^2-x^2),  color=254, linestyle=1
     OPLOT, x, -SQRT((i*sol_rad/xscale)^2-x^2), color=254, linestyle=1
   ENDFOR
   FOR i=10,25, 5 DO BEGIN
     x=FINDGEN(arr)*(i*sol_rad/xscale)/DOUBLE(arr-1)
     OPLOT, x, SQRT((i*sol_rad/xscale)^2-x^2),  color=254, linestyle=1
     OPLOT, x, -SQRT((i*sol_rad/xscale)^2-x^2), color=254, linestyle=1
   ENDFOR
   FOR i=50,200, 50 DO BEGIN
     x=FINDGEN(arr)*(i*sol_rad/xscale)/DOUBLE(arr-1)
     OPLOT, x, SQRT((i*sol_rad/xscale)^2-x^2),  color=254, linestyle=1
     OPLOT, x, -SQRT((i*sol_rad/xscale)^2-x^2), color=254, linestyle=1
   ENDFOR

   ;______________
   ;
   ; Add colourbar 
   ;______________ 

   COLORBAR,  CHARSIZE = 1.5, DIVISIONS = 5, FORMAT = '(F4.2)',  $ 
   POSITION = [0.23, 0.023, 0.38, 0.043],  $ 
   NCOLORS = NGLev-1, RANGE = [minkappa, maxkappa], BOTTOM = 2
  ;--------------------------------------------------------------------

   IF(nos GT 0)THEN BEGIN
     NGLev=nos
  ;--------------------------------------------------------------------
    ;get a previously saved colour-map
     restore, 'colours/Structure'+nos_string+'.sav' & tvlct, r, g, b            

     ;set levels for contour plot
     stepg=nos/DOUBLE(NGLev) & grey=indgen(NGLev)*stepg+1;

     ;plot the structures
     contour, Structures[*,*,0], (FIndGen(XHs)*RangeGm)/DOUBLE(XHs-1),    $ 
              ((FIndGen(2*YHs-1)*RangeGm)/DOUBLE(YHs-1)-RangeGm),         $ 
              TITLE = TITLE_FONT+'structures !3',                         $ 
              XTitle='X!DH!N  (Gm)', YTitle='Y!DH!N  (Gm)', charsize=1.7, $ 
              levels=grey, /fill, c_colors=1+indgen(NGLev)
     FOR i=1,nos-1 DO BEGIN
       contour, Structures[*,*,i], (FIndGen(XHs)*RangeGm)/DOUBLE(XHs-1),  $ 
                ((FIndGen(2*YHs-1)*RangeGm)/DOUBLE(YHs-1)-RangeGm),       $
                levels=grey, /fill, c_colors=1+indgen(NGLev),             $
               /OVERPLOT
     ENDFOR ;plot the structures

    ;Sun
     x=FINDGEN(arr)*(sol_rad/xscale)/DOUBLE(arr-1)
     OPLOT, x, SQRT((sol_rad/xscale)^2-x^2),  color=0, THICK=4
     OPLOT, x, -SQRT((sol_rad/xscale)^2-x^2), color=0, THICK=4 

     ;--------------------------------------------------------------------

      ;______________
      ;
      ; Add colourbar 
      ;______________ 

      COLORBAR, CHARSIZE = 1.5,             $ 
                POSITION = [0.43, 0.023, 0.58, 0.043],            $ 
                FORMAT = '(I1)', DIVISIONS = nos, RANGE = [0, nos], $ 
                NCOLORS = NGLev, BOTTOM = 1
     ;--------------------------------------------------------------------     
   ENDIF   ;(nos GT 0)
   ;-------------------------------------------------------------------- 

   NCLev=7 
   NGLev=8*NCLev  
   restore, 'colours/DS_WG.sav'  & tvlct, r, g, b 

  ;set levels for contour plot
   stepg=HF_Range/DOUBLE(NGLev);
   grey=indgen(NGLev)*stepg+MinHF;

  ;Fundamental Emission intensities
   contour, Fund, (FIndGen(XH)*Radial_Height[t_steps-1])/DOUBLE(XH-1)/xscale,   $
                  ((FIndGen(2*YH-1)*Radial_Height[t_steps-1])/DOUBLE(YH-1)      $
                   -Radial_Height[t_steps-1])/yscale,                           $
                  Title=TITLE_FONT+' Fundamental !3',                           $
                  XTitle='X!DH!N  (Gm)', YTitle='Y!DH!N  (Gm)', charsize=1.7,   $
                   XRange=[0,RangeGm], YRange=[-RangeGm,RangeGm],               $
                   levels=grey, /fill, c_colors=1+indgen(NGLev)

  ;Observer Location  
   xyouts, observerX, observerY-RangeGm/15, '!4!I*!N!3',  $
           color=55, charsize=7.0, ALIGNMENT=0.5

  ;Global shocks at periodic time steps
   FOR i=0,t_steps-1, k DO BEGIN
   IF((Radial_Height[i] GT sol_rad) AND (Radial_Height[i] LT RangeGm*1.0e9))THEN BEGIN
    x_shk=FINDGEN(arr)*(Radial_Height[i])/DOUBLE(arr-1)
    Y_shk=SQRT(x_shk/CMEb[i]) & Y_shk= ROTATE(Y_shk,2)
    Y2_shk= -Y_shk
    j=LONG(sol_rad/(Radial_Height[i]/DOUBLE(arr-1)))
    Plots, x_shk[j:(arr-1)]/xscale, Y_shk[j:(arr-1)]/yscale, THICK=3 
    Plots, x_shk[j:(arr-1)]/xscale, Y2_shk[j:(arr-1)]/yscale, THICK=3
   ENDIF
   ENDFOR

  ;Sun
   x=FINDGEN(arr)*(sol_rad/xscale)/DOUBLE(arr-1)
   OPLOT, x, SQRT((sol_rad/xscale)^2-x^2),  color=254, THICK=4
   OPLOT, x, -SQRT((sol_rad/xscale)^2-x^2), color=254, THICK=4

  ;Concentric circles at integer numbers of solar radii
   FOR i=2,5 DO BEGIN
     x=FINDGEN(arr)*(i*sol_rad/xscale)/DOUBLE(arr-1)
     OPLOT, x, SQRT((i*sol_rad/xscale)^2-x^2),  color=254, linestyle=1 
     OPLOT, x, -SQRT((i*sol_rad/xscale)^2-x^2), color=254, linestyle=1 
   ENDFOR
   FOR i=10,25, 5 DO BEGIN
     x=FINDGEN(arr)*(i*sol_rad/xscale)/DOUBLE(arr-1)
     OPLOT, x, SQRT((i*sol_rad/xscale)^2-x^2),  color=254, linestyle=1 
     OPLOT, x, -SQRT((i*sol_rad/xscale)^2-x^2), color=254, linestyle=1 
   ENDFOR
   FOR i=50,200, 50 DO BEGIN
     x=FINDGEN(arr)*(i*sol_rad/xscale)/DOUBLE(arr-1)
     OPLOT, x, SQRT((i*sol_rad/xscale)^2-x^2),  color=254, linestyle=1 
     OPLOT, x, -SQRT((i*sol_rad/xscale)^2-x^2), color=254, linestyle=1 
   ENDFOR

  ;--------------------------------------------------------------------

  ;Harmonic Emission intensities
   contour, Harm, (FIndGen(XH)*Radial_Height[t_steps-1])/DOUBLE(XH-1)/xscale,  $ 
                  ((FIndGen(2*YH-1)*Radial_Height[t_steps-1])/DOUBLE(YH-1)     $ 
                  -Radial_Height[t_steps-1])/yscale,                           $
                  Title=TITLE_FONT+' Harmonic !3',                             $
                  XTitle='X!DH!N  (Gm)', YTitle='Y!DH!N  (Gm)', charsize=1.7,  $ 
                  XRange=[0,RangeGm], YRange=[-RangeGm,RangeGm],               $
                  levels=grey, /fill, c_colors=1+indgen(NGLev)

  ;Observer Location  
   xyouts, observerX, observerY-RangeGm/15, '!4!I*!N!3',  $ 
           color=55, charsize=7.0, ALIGNMENT=0.5

  ;Global shocks at periodic time steps
   FOR i=0,t_steps-1, k DO BEGIN
   IF((Radial_Height[i] GT sol_rad) AND (Radial_Height[i] LT RangeGm*1.0e9))THEN BEGIN
    x_shk=FINDGEN(arr)*Radial_Height[i]/DOUBLE(arr-1)
    Y_shk=SQRT(x_shk/CMEb[i]) & Y_shk= ROTATE(Y_shk,2)
    Y2_shk= -Y_shk
    j=LONG(sol_rad/(Radial_Height[i]/DOUBLE(arr-1)))   
    Plots, x_shk[j:(arr-1)]/xscale, Y_shk[j:(arr-1)]/yscale, THICK=3 
    Plots, x_shk[j:(arr-1)]/xscale, Y2_shk[j:(arr-1)]/yscale, THICK=3
   ENDIF
   ENDFOR

  ;Sun
   x=FINDGEN(arr)*(sol_rad/xscale)/DOUBLE(arr-1)
   OPLOT, x, SQRT((sol_rad/xscale)^2-x^2),  color=254, THICK=4
   OPLOT, x, -SQRT((sol_rad/xscale)^2-x^2), color=254, THICK=4

  ;Concentric circles at integer numbers of solar radii
   FOR i=2,5 DO BEGIN
     x=FINDGEN(arr)*(i*sol_rad/xscale)/DOUBLE(arr-1)
     OPLOT, x, SQRT((i*sol_rad/xscale)^2-x^2),  color=254, linestyle=1 
     OPLOT, x, -SQRT((i*sol_rad/xscale)^2-x^2), color=254, linestyle=1 
   ENDFOR
   FOR i=10,25, 5 DO BEGIN
     x=FINDGEN(arr)*(i*sol_rad/xscale)/DOUBLE(arr-1)
     OPLOT, x, SQRT((i*sol_rad/xscale)^2-x^2),  color=254, linestyle=1 
     OPLOT, x, -SQRT((i*sol_rad/xscale)^2-x^2), color=254, linestyle=1 
   ENDFOR
   FOR i=50,200, 50 DO BEGIN
     x=FINDGEN(arr)*(i*sol_rad/xscale)/DOUBLE(arr-1)
     OPLOT, x, SQRT((i*sol_rad/xscale)^2-x^2),  color=254, linestyle=1 
     OPLOT, x, -SQRT((i*sol_rad/xscale)^2-x^2), color=254, linestyle=1 
   ENDFOR

  ;--------------------------------------------------------------------

   ;______________
   ;
   ; Add colourbar 
   ;______________ 

   COLORBAR, TITLE = 'Flux (W m!U-2!N Sr!U-1!N)', CHARSIZE = 1.5,            $ 
   POSITION = [0.63, 0.023, 0.98, 0.043], FORMAT = '(F5.1)', DIVISIONS = 5,  $ 
   NCOLORS = NGLev-1, RANGE = [minHF, maxHF], BOTTOM = 2

  ;__________________________________________________________________

   CloseDevice, OPFilename, OPDevice
 endif                
;--------------------------------------------------------------------     

end ;Helio_cont.pro

; consider: bin all calculated ripples from all time steps as a function 
;           of radial height then create the 3D data by filling the 
;           annulus of each ripple with randomly chosen examples from 
;           its radial heigth

