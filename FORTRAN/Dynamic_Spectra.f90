!Dynamic_Spectra.f90  (Thesis Version)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Fortran 90    First Written: 18/12/2001 (Reorganization/extension of 
!!                                           a code begun in 1999)   
!!                Current Version: see MODIFICATION HISTORY below
!!  
!!  Written by: Stuart Knock 
!!
!!  PURPOSE: Calculates the predicted dynamic spectra for a type II  
!!           radio burst from the solar corona and interplanetary 
!!           medium.
!!
!!  MAIN PROGRAM: Dynamic_Spectra --> Loops over the ripples which make 
!!                                    up a global shock front for 
!!                                    successive time steps, calling 
!!                                    relevant subroutines while deciding 
!!                                    whether to continue and writing 
!!                                    some information to file...
!!                                    
!!
!!  MODULES: Global --> kind definitions and physical parameters used 
!!                      throughout the code.
!!
!!           Shock --> parameter definitions and declarations related to 
!!                     the shock. eg curvature, strength, speed etc.
!!
!!           Wind --> parameter definitions and declarations related to
!!                    the solar wind. eg speed, temperature, density etc.
!!
!!           Beam --> parameter definitions and declarations related to 
!!                    the electron beams. eg number density, width etc.
!!                    CONTAINS gamma function FUNCTION
!!
!!           Emis --> parameter definitions and declarations related to 
!!                    the calculation of volume emissivities. eg solid 
!!                    angle over which emission is spread.
!!                    CONTAINS error function FUNCTION
!!
!!           OutPut --> file name and other variables associated with the 
!!                      programs output.  
!!
!!           Spectra --> Arrays & associated variables of the calculated 
!!                       volume emissivities. 
!!
!!           NR --> contains Numerical Recipes routines including: 
!!                  - Gaussian_Random(harvest)
!!
!!           
!!  SUBROUTINES: Initialise_Global --> Reads in initial parameters and 
!!                                     ranges from Dynamic_Spectra.input
!!                                     and writes some information about 
!!                                     initial values to a LOG_*.dat file.
!!                                     CALLED FROM MAIN.
!!               
!!               Initialise_Time_Step --> sets up things that need to be 
!!                                        set at the beginning of each
!!                                        time step. CALLED FROM MAIN.
!!                                        
!!               
!!               Ripple_Parameters --> Determines the shock and upstream 
!!                                     conditions for the ripples that 
!!                                     make up the global shock. CALLED 
!!                                     FROM MAIN.
!!               
!!               bckgrnd_distribution --> calculates the reduced solar 
!!                                        wind (background) distribution 
!!                                        function for a given kappa and
!!                                        electron thermal speed Ve.
!!                                        CALLED FROM RIPPLE_PARAMETERS.
!!               
!!               shock_parameter_arrays --> calculates shock parameters 
!!                                          as a function of location on 
!!                                          a shock ripple. CALLED FROM 
!!                                          RIPPLE_PARAMETERS.
!!               
!!               B2B1_ephi --> calculates the magnetic field jump and 
!!                             cross-shock-potential  as a function of 
!!                             position on a ripple. CALLED FROM 
!!                             SHOCK_PARAMETER_ARRAYS. 
!!               
!!               cgroots --> modified numerical recipes routine to set 
!!                           up and call ZROOTS & LAGUER to determine 
!!                           density jump at a shock from 
!!                           Rankine-Hugoniot conditions. CALLED FROM 
!!                           B2B1_EPHI.
!!               
!!               ZROOTS --> density jump. CALLED FROM CGROOTS.
!!               
!!               LAGUER --> density jump. CALLED FROM ZROOTS.
!!               
!!               Flux --> Calls the routines that calculate the emission 
!!                        from a given ripple and assigns the emission to 
!!                        dynamic spectra arrays for output. CALLED
!!                        FROM MAIN.
!!               
!!               cutoff_velocity_distribution --> Calculates the 
!!                                                reflected component of 
!!                                                the foreshock 
!!                                                distribution function 
!!                                                for a given foreshock 
!!                                                location. CALLED FROM 
!!                                                FLUX.
!!               
!!               qlf --> Takes a beam distribution, created by combining 
!!                       a cutoff distribution with the background wind
!!                       distribution, and determines the beam parameters 
!!                       (ie width, speed, density). CALLED FROM FLUX.
!!               
!!               emissivities --> Takes in beam parameters an calculates 
!!                                fundamental and harmonic volume 
!!                                emissivities. CALLED FROM FLUX.
!!               
!!               Finalise_and_Write --> Sets up dynamic spectra for 
!!                                      output to a file, and then 
!!                                      clears the output data from 
!!                                      stored arrays. CALLED FROM MAIN.
!!               
!!               Write_to_file --> CALLED FROM FINALISE_AND_WRITE. 
!!               
!!               
!!  USAGE: see README file in this directory
!!               
!!  INPUTS: see Dynamic_Spectra.input & README file in this directory
!!
!!  OUTPUTS: see README file in this directory                      
!!                  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MODIFICATION HISTORY (Most recent modification first)
!! date is when modification is basically functioning, 
!! though further refinement may be made afterward.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!echo -n \!! \ Last Significant Modification On \ && date && echo -n \!! \ 
!!  Last Significant Modification On  Thu Feb 19 12:38:36 EST 2004
!!  
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!____________________________________________________________________
! 
!                             MODULES 
!____________________________________________________________________
 
  MODULE Global
! *************
  IMPLICIT NONE
  SAVE
!
! Define Precision:
!
  INTEGER, PARAMETER :: L=SELECTED_REAL_KIND(13,300) !(sig digs,range)
!
! Declare physical constants:
! c=speed of light (m/s)            mi=ion mass (assuming protons)(kg)
! kbolt=Boltzman constant           me=electron mass (kg)
! AU=Sun-Earth distance (m)         Sol_rad=Radius of Sun (m)
! Sol_period=Solar period (s)
! 
  REAL(KIND=L), PARAMETER :: c=2.9979e8_L, kbolt=1.380622e-23_L,   &
                             me=9.109558e-31_L, mi=1.672614e-27_L, &
                             AU=1.49599e11_L, Sol_rad=6.965e8_L,   & 
                             Sol_period=2194560.0_L 
!
  REAL(KIND=L), PARAMETER :: me2on=2.0_L/me   
!  
! Other constants used frequently:  
!  
  REAL(KIND=L), PARAMETER :: pi=3.141592653589793238462643_L 
  REAL(KIND=L), PARAMETER :: mu0=4.0_L*pi*1.0e-7_L
  REAL(KIND=L), PARAMETER :: pi15=1.5_L*pi 
  REAL(KIND=L), PARAMETER :: rtpi=1.77245385090552_L
  REAL(KIND=L), PARAMETER :: sqthree=1.73205080756888_L
  REAL(KIND=L), PARAMETER :: sqtwo=1.41421356237310_L
  REAL(KIND=L), PARAMETER :: sqpis=0.723601254558268_L
  REAL(KIND=L), PARAMETER :: csq=c*c
!
  INTEGER :: calculate
!
  END MODULE Global
!____________________________________________________________________

  MODULE Shock
! ************
  USE Global, ONLY: L
  IMPLICIT NONE
  SAVE
!
  REAL(KIND=L), PARAMETER :: Avg_num_per_timestep=100.0_L
  INTEGER, PARAMETER :: SR=30000           ! shock param. arr. size
  INTEGER, PARAMETER :: RR=100, xR=200     ! dimensions wrt B-field
!
  REAL(KIND=L), DIMENSION(SR) :: alph=0.0_L    ! Ang. surf.->Bfield
  REAL(KIND=L), DIMENSION(SR) :: thetabn=0.0_L ! Ang. Bfield->norm.
  REAL(KIND=L), DIMENSION(SR) :: Vc=0.0_L      ! cutoff vel 
  REAL(KIND=L), DIMENSION(SR) :: Un=0.0_L      ! normal component U
  REAL(KIND=L), DIMENSION(SR) :: UnSalph=0.0_L ! Un*SIN(alph)
  REAL(KIND=L), DIMENSION(SR) :: n2n1=0.0_L    ! numb. dens. ratio
  REAL(KIND=L), DIMENSION(SR) :: Ma=0.0_L      ! Alfven Mach number
  REAL(KIND=L), DIMENSION(SR) :: Mssq=0.0_L    ! Sonic-Mach-number^2
  REAL(KIND=L), DIMENSION(SR) :: B2B1=0.0_L    ! B-field ratio 
  REAL(KIND=L), DIMENSION(SR) :: B2B1_1=0.0_L  ! B-field ratio -1 
  REAL(KIND=L), DIMENSION(SR) :: UnSalphonB2B1=0.0_L
  REAL(KIND=L), DIMENSION(SR) :: ephi=0.0_L    ! cross-shock poten.
  REAL(KIND=L), DIMENSION(SR) :: ephi2onme=0.0_L
!
  REAL(KIND=L) :: Radial_Height  ! height; leading point in meters
  REAL(KIND=L) :: Shock_Stop     ! height where calculation stops
  REAL(KIND=L) :: Ripple_Height  ! height of ripple above sun in meters
  REAL(KIND=L) :: RHonAU         ! Ripple height as fraction of AU
  REAL(KIND=L) :: CME_b          ! constant of curvature; global
  REAL(KIND=L) :: b              ! constant of curvature; ripple 
  REAL(KIND=L) :: Average_b    
  REAL(KIND=L) :: U              ! speed wrt wind
  REAL(KIND=L) :: R, x           ! coordinates parallel & perp to B 
  REAL(KIND=L) :: Rrange, xrange ! extent of source region calculation
  REAL(KIND=L) :: Srange         ! range of shock R in calc region
  REAL(KIND=L) :: theta_radial, alpha_global 
  REAL(KIND=L) :: RHSX
! 
  INTEGER :: UT_time    ! current time in simulation 
  INTEGER :: t_step     ! time resolution in seconds 
  INTEGER :: S_limit
!  
! Intermediate calculation variables:
!
  REAL(KIND=L) :: SRonSrange
  REAL(KIND=L) :: UCME, CME_accel  ! Shock Speed wrt Sun, deceleration
  REAL(KIND=L) :: rip_pck_num, rip_spacing, SX, SY, SY_temp
! 
  REAL(KIND=L) :: expansion 
  REAL(KIND=L) :: scale, scale_exp
! 
  END MODULE Shock
!____________________________________________________________________

  MODULE Wind
! ***********
  USE Global, ONLY: L
  IMPLICIT NONE
  SAVE
!
! Structure definition for structures in the solar wind:
! 
  TYPE SW_Struc
    INTEGER :: shape
    INTEGER :: nop
    REAL(KIND=L), ALLOCATABLE, DIMENSION(:) :: param
    REAL(KIND=L) :: theta,Ne,B1,Ti,Te,Vsw,b,Kappa
  END TYPE SW_Struc 
!
  TYPE(SW_Struc), ALLOCATABLE, DIMENSION(:) :: SWS
!
  INTEGER :: nos  !number of structures
!
! Declare physical variables: 
!
  REAL(KIND=L) :: Te,Ti       ! electron & ion temperatures
  REAL(KIND=L) :: Average_Te,Average_Ti
  REAL(KIND=L) :: Te_1AU,Ti_1AU      
  REAL(KIND=L) :: Ne          ! electron number density
  REAL(KIND=L) :: Average_Ne 
  REAL(KIND=L) :: Average_fp 
  REAL(KIND=L) :: Ne_1AU,Ne_1Rs          
  REAL(KIND=L) :: B1          ! Upstream magnetic field strength
  REAL(KIND=L) :: Average_B1      
  REAL(KIND=L) :: Bo         
  REAL(KIND=L) :: Vsw         ! Solar Wind Speed
  REAL(KIND=L) :: Average_Vsw        
  REAL(KIND=L) :: Vsw_1AU  
  REAL(KIND=L) :: theta_Parker
  REAL(KIND=L) :: theta       !  Ang. -U -> Bfield  in radians
  REAL(KIND=L) :: Average_theta
  REAL(KIND=L) :: Vd          ! drift velocity (E X B)
  REAL(KIND=L) :: Ve          ! electron thermal velocity
  REAL(KIND=L) :: cs          ! sound speed
  REAL(KIND=L) :: kappa       ! e- distribution function parameter
  REAL(KIND=L) :: Average_kappa 
!
! Intermediate calculation variables:
!
  REAL(KIND=L) :: stheta, stheta2, ctheta, tantheta, stheta1on, Vesq
!
  END MODULE Wind
!____________________________________________________________________

  MODULE Beam
! ***********
  USE Global, ONLY: L
  IMPLICIT NONE
  SAVE
!
  INTEGER, PARAMETER :: VR=20000    ! beam array size
!  
  REAL(KIND=L), DIMENSION(VR) :: vparray=0.0_L ! vel|| dist'bn
  REAL(KIND=L), DIMENSION(VR) :: vcarray=0.0_L ! cutoff vel dist'bn
  REAL(KIND=L), DIMENSION(VR) :: f=0.0_L         ! total vel dist'bn
  REAL(KIND=L), DIMENSION(VR) :: vel=0.0_L     ! beam velocity
  REAL(KIND=L), DIMENSION(VR) :: temp4vponvd=0.0_L
!  
! Maximum velocity to which distributions are calculated:
  REAL(KIND=L), PARAMETER :: vmax=2.0e8_L
  REAL(KIND=L), PARAMETER :: vmaxonVR=vmax/VR
!  
! beam velocity, beam width  and number of electrons in beam:
  REAL(KIND=L) :: vb, deltavb, Nb
  REAL(KIND=L) :: dv, halfdv
!
! Intermediate calculation variables:
!
  INTEGER :: vcnt, vend, vendold, vcntold 
  REAL(KIND=L) :: ConstDist,  b4on, bsqsthetasq, bstheta1on, Uctp
! 
!
  CONTAINS 
   REAL(KIND=L) FUNCTION gammafn(xx)
!  **************************************
   IMPLICIT NONE
!
!   Calculates a gamma function value for xx
!    (a hack from numerical recipies)
!
    INTEGER, PARAMETER :: dp = KIND(1.0D0) 
    REAL(KIND=L), INTENT(IN) :: xx
    REAL(KIND=dp) :: lngammafn
    REAL(KIND=dp) :: tmp, xxx
    REAL(KIND=dp), PARAMETER :: stp = 2.5066282746310005_dp 
    REAL(KIND=dp), DIMENSION(6) :: arthp
    REAL(KIND=dp), DIMENSION(6) :: coef = (/76.18009172947146_dp, &
                     -86.50532032941677_dp, 24.01409824083091_dp, & 
                  -1.231739572450155_dp,0.1208650973866179e-2_dp, & 
                  -0.5395239384953e-5_dp/)
    INTEGER :: ii
!  
    xxx=xx 
    tmp=xxx+5.5_dp
    tmp=(xxx+0.5_dp)*log(tmp)-tmp
    arthp(1)=xxx+1.0_dp
    DO ii=2,6 ! 6 is SIZE(coef)
      arthp(ii)=arthp(ii-1)+1.0_dp
    END DO
    lngammafn=tmp+log(stp*(1.000000000190015_dp+SUM(coef/arthp))/xxx)
    gammafn=EXP(lngammafn)
!
   END FUNCTION gammafn
!
  END MODULE Beam
!____________________________________________________________________

  MODULE Emis
! ***********
  USE Global, ONLY: L, pi 
  IMPLICIT NONE
  SAVE
!  
! Solid angle spread of fundamental and harmonic radio emission:
! 
  REAL(KIND=L), PARAMETER :: delohmF=0.25_L*pi, delohmH=2.0_L*pi 
! 
! Intermediate calculation variables:
!
  REAL(KIND=L) :: temp4gammaLS, Vcrit, temp4zetaF, erftemp1, & 
                  erftemp2, temp4phiF, temp4phiH, gmemisq 
!
  CONTAINS 
   REAL(KIND=L) FUNCTION erf(x)
!  *******************************
   IMPLICIT NONE
!
!  Calculates error function on real axis.
!  Accuracy = 1.5e-7.
!  Reference: Abramowitz M. and Stegun I. A. p.299, eq. 7.1.26
!
   REAL(KIND=L),INTENT(IN) :: x 
   REAL(KIND=L) :: t,tsq
   REAL(KIND=L),PARAMETER :: a1=0.254829592_L,a2=-0.284496736_L, &
          a3=1.421413741_L,a4=-1.453152027_L,a5=1.061405429_L,   & 
          p=0.3275911_L
!
!  print*,'function being evaluated=',x
!
   IF(x > 99.0_L)THEN
     erf=1.0_L
   ELSE IF(x < -99.0_L)THEN
     erf=-1.0_L
   ELSE  
     t=1.0_L/(1.0_L+p*ABS(x))
     tsq=t*t
     erf=a1*t+a2*tsq+a3*t*tsq+a4*tsq*tsq+a5*tsq*tsq*t
     erf=1.0_L-erf*exp(-x*x)
     IF(x < 0.0_L)erf=-erf
   END IF    
!
   END FUNCTION erf
!
   SUBROUTINE Gaussian_Random_e(harvest)
!  ***********************************
   USE Global, ONLY: L
   IMPLICIT NONE 
   REAL(KIND=L), DIMENSION(:), INTENT(OUT) :: harvest
   REAL(KIND=L), DIMENSION(SIZE(harvest)) :: rsq, v1, v2 
   REAL(KIND=L), ALLOCATABLE, DIMENSION(:), SAVE :: g 
   INTEGER :: n, i
   LOGICAL, SAVE :: gaus_stored=.FALSE. 
!
   n=SIZE(harvest)
   IF(.NOT. gaus_stored)ALLOCATE(g(n))   
!  
   IF(gaus_stored)THEN     
     harvest=g            
     gaus_stored=.FALSE.  
   ELSE                   
     DO i=1,n
       DO 
         CALL RANDOM_NUMBER(v1(i))
         CALL RANDOM_NUMBER(v2(i))
         v1(i)=2.0_L*v1(i)-1.0_L
         v2(i)=2.0_L*v2(i)-1.0_L
         rsq(i)=v1(i)*v1(i)+v2(i)*v2(i) 
         IF(rsq(i) > 0.0 .AND. rsq(i) < 1.0)EXIT
       END DO
     END DO
     rsq=SQRT(-2.0_L*LOG(rsq)/rsq)
     harvest=(v1*rsq)*0.25_L  !divide by 4 makes  SD 0.5
     g=(v2*rsq)*0.25_L        ! or by 6 makes  SD 0.333...
     gaus_stored=.TRUE.
   END IF
!
   IF(.NOT. gaus_stored)DEALLOCATE(g) 
!
   END SUBROUTINE Gaussian_Random_e  
!
  END MODULE Emis
!____________________________________________________________________ 
  
  MODULE OutPut
! **************
  USE Global, ONLY: L
  IMPLICIT NONE
  SAVE
! 
  INTEGER, PARAMETER  :: PR=100 !Array size of parameters
  INTEGER, DIMENSION(PR) :: theta_arr, U_arr, b_arr, B1_arr, & 
                            Ne_arr, Te_arr, Ti_arr
  REAL(KIND=L) :: total_H1, total_F1
  REAL(KIND=L) :: total_H2, total_F2 
  REAL(KIND=L) :: total_H3, total_F3
!
! Parameter plot ranges
!
  REAL(KIND=L) :: Maxtheta, Mintheta
  REAL(KIND=L) :: MaxU, MinU
  REAL(KIND=L) :: Maxb, Minb
  REAL(KIND=L) :: MaxB1, MinB1
  REAL(KIND=L) :: MaxNe, MinNe
  REAL(KIND=L) :: MaxTe, MinTe
  REAL(KIND=L) :: MaxTi, MinTi
!
  INTEGER :: first_pass !used to signify 1st time writing to file
  CHARACTER(LEN=43) :: Output_dir
  CHARACTER(LEN=17) :: file_ident !common ending for various output files
!                           
  END MODULE OutPut              
!____________________________________________________________________
   
  MODULE Spectra                                          
! **************                                                
  USE Global, ONLY: L
  IMPLICIT NONE
  SAVE
! 
  INTEGER, PARAMETER :: freq_bins=1000    ! bins in frequency
  INTEGER, PARAMETER :: time_bins=1440    ! bins in time_window
  INTEGER, PARAMETER :: time_window=14400 ! time_window in seconds
  INTEGER, PARAMETER :: time_resolution=time_window/time_bins  !of DS in seconds 
  REAL(KIND=L), DIMENSION(freq_bins,time_bins) :: DS1F, DS2F, DS3F,  & 
                                                  DS1H, DS2H, DS3H
  REAL(KIND=L) :: exp_min, exp_range, l_val
  REAL(KIND=L) :: fber, erfb  ! exponent range to frequency bin ratios
  REAL(KIND=L) :: loglem
  REAL(KIND=L) :: Xo1,Xo2,Xo3
  REAL(KIND=L) :: Yo1,Yo2,Yo3
  INTEGER :: spec_shift, blocked1,blocked2,blocked3
!
  END MODULE Spectra
!____________________________________________________________________
   
  MODULE NR                                          
! *********                                                
  USE Global, ONLY: L
  IMPLICIT NONE
  SAVE  
! 
  CONTAINS 
   SUBROUTINE Gaussian_Random(harvest)
!  ***********************************
   USE Global, ONLY: L
   IMPLICIT NONE 
   REAL(KIND=L), DIMENSION(:), INTENT(OUT) :: harvest
   REAL(KIND=L), DIMENSION(SIZE(harvest)) :: rsq, v1, v2 
   REAL(KIND=L), ALLOCATABLE, DIMENSION(:), SAVE :: g 
   INTEGER :: n, i
   LOGICAL, SAVE :: gaus_stored=.FALSE. 
!
   n=SIZE(harvest)
   IF(.NOT. gaus_stored)ALLOCATE(g(n))   
!  
   IF(gaus_stored)THEN    
     harvest=g            
     gaus_stored=.FALSE.  
   ELSE 
     DO i=1,n
       DO 
         CALL RANDOM_NUMBER(v1(i))
         CALL RANDOM_NUMBER(v2(i))
         v1(i)=2.0_L*v1(i)-1.0_L
         v2(i)=2.0_L*v2(i)-1.0_L
         rsq(i)=v1(i)*v1(i)+v2(i)*v2(i) 
         IF(rsq(i) > 0.0 .AND. rsq(i) < 1.0)EXIT
       END DO
     END DO
     rsq=SQRT(-2.0_L*LOG(rsq)/rsq)
     harvest=(v1*rsq)*0.25_L  !divide by 4 makes  Standard Deviation 0.5
     g=(v2*rsq)*0.25_L        ! or by 6 makes  SD 0.333...
     gaus_stored=.TRUE.
   END IF
!
   IF(.NOT. gaus_stored)DEALLOCATE(g) 
!
   END SUBROUTINE Gaussian_Random  
!
  END MODULE NR
!____________________________________________________________________

!____________________________________________________________________
! 
!                           MAIN PROGRAM 
!____________________________________________________________________
 
!********************************************************************
  PROGRAM Dynamic_Spectra
!******************************************************************** 
  USE Global, ONLY: L, AU, Sol_rad, pi, calculate
  USE Shock, ONLY: CME_b, Radial_Height, Shock_Stop, alpha_global,    & 
                   SX, SY, RHSX, UCME, Avg_num_per_timestep, SY_temp, & 
                   rip_pck_num, rip_spacing, b, Average_b
  USE OutPut, ONLY: total_F1, total_H1, total_F2, total_H2, & 
                    total_F3, total_H3, file_ident, Output_dir
  USE wind, ONLY: Average_fp, theta, Ne, B1, Ti, Te, Vsw, Kappa
  IMPLICIT NONE
  INTEGER :: i, ripples 
  CHARACTER(LEN=8) :: DATE
  CHARACTER(LEN=10) :: TIME
666 FORMAT('END DATE: ',A2,'/',A2,'/',A4,'     ',  & 
           'END TIME: ',A2,':',A2,':',A2,/)
! 
  CALL Initialise_Global
!
  OPEN(unit=23,file=TRIM(Output_dir)//'Source_Loc1_'//file_ident,status='new')
  OPEN(unit=24,file=TRIM(Output_dir)//'Source_Loc2_'//file_ident,status='new')
  OPEN(unit=25,file=TRIM(Output_dir)//'Source_Loc3_'//file_ident,status='new')
  OPEN(unit=37,file=TRIM(Output_dir)//'Wind_Loc_'//file_ident,status='new')
  OPEN(unit=27,file=TRIM(Output_dir)//'Source_Num_'//file_ident,status='new') 
  ENDFILE(unit=27) 
  CLOSE(unit=27)
  OPEN(unit=29,file=TRIM(Output_dir)//'CME_b_'//file_ident,status='new')
  ENDFILE(unit=29) 
  CLOSE(unit=29)
!
  DO !time_step            
    CALL Initialise_Time_Step
    ripples=0
    OPEN(unit=29, file=TRIM(Output_dir)//'CME_b_'//file_ident, & 
         status='old', position='append')
      write(unit=29,fmt=*)CME_b
    ENDFILE(unit=29) 
    CLOSE(unit=29)
    IF(Radial_Height >= Shock_Stop)THEN 
      OPEN(unit=21, file=TRIM(Output_dir)//'LOG_'//file_ident, & 
           status='old', position='append') 
        CALL DATE_AND_TIME(DATE,TIME) 
        write(unit=21,fmt=666)DATE(7:8),DATE(5:6),DATE(1:4), & 
                              TIME(1:2),TIME(3:4),TIME(5:6)
        write(unit=21,fmt=*)'Exceeded Radial_Height limit'
        ENDFILE(unit=21) 
      CLOSE(unit=21) 
      EXIT
    END IF  
    IF(UCME <= 300.0e3)THEN 
      OPEN(unit=21, file=TRIM(Output_dir)//'LOG_'//file_ident, & 
           status='old', position='append') 
        CALL DATE_AND_TIME(DATE,TIME)
        write(unit=21,fmt=666)DATE(7:8),DATE(5:6),DATE(1:4), & 
                              TIME(1:2),TIME(3:4),TIME(5:6) 
        write(unit=21,fmt=*)'Low UCME'
        ENDFILE(unit=21) 
      CLOSE(unit=21) 
      EXIT 
    END IF
    DO !ripples on global shock
      DO i=1,2 !either side
        total_F1=0.0_L ; total_H1=0.0_L
        total_F2=0.0_L ; total_H2=0.0_L
        total_F3=0.0_L ; total_H3=0.0_L
        CALL Ripple_Parameters (i)
        write(unit=37,fmt=*)RHSX,SY_temp,theta,Ne,B1,Ti,Te,Vsw,Kappa
        IF(calculate == 1)CALL Flux 
        calculate=1
        ripples=ripples+1
        write(unit=23,fmt=*)RHSX,SY_temp,total_F1,total_H1,Average_fp
        write(unit=24,fmt=*)RHSX,SY_temp,total_F2,total_H2,Average_fp
        write(unit=25,fmt=*)RHSX,SY_temp,total_F3,total_H3,Average_fp
      END DO !either side
      rip_spacing=3.0_L/b
      SY=SY+COS(alpha_global)*rip_spacing
      SX=CME_b*SY*SY
      rip_pck_num=SY*pi/rip_spacing !semi circle of global 
      IF(SX >= 0.95*(Radial_Height-Sol_rad))EXIT
    END DO !ripples on global shock
    CALL Finalise_and_Write
    OPEN(unit=27, file=TRIM(Output_dir)//'Source_Num_'//file_ident, & 
         status='old', position='append')
      write(unit=27,fmt=*)ripples
    ENDFILE(unit=27) 
    CLOSE(unit=27)
  END DO !time_step
!
  ENDFILE(unit=23) 
  CLOSE(unit=23)
  ENDFILE(unit=24) 
  CLOSE(unit=24)
  ENDFILE(unit=25) 
  CLOSE(unit=25)
  ENDFILE(unit=37) 
  CLOSE(unit=37)
!
  END PROGRAM Dynamic_Spectra
!********************************************************************
!____________________________________________________________________
!

!____________________________________________________________________
! 
!                           SUBROUTINES 
!____________________________________________________________________

  SUBROUTINE Initialise_Global
! ****************************    
  USE Global, ONLY: L, Sol_rad
  USE OutPut, ONLY: first_pass, file_ident, Output_dir, PR,            & 
                    Maxtheta, Mintheta, MaxU, MinU, Maxb, Minb, MaxB1, & 
                    MinB1, MaxNe, MinNe, MaxTe, MinTe, MaxTi, MinTi
  USE Shock, ONLY: UT_time, t_step, Radial_Height, Shock_Stop, UCME, & 
                   CME_accel, SR, RR, xR, Avg_num_per_timestep,      & 
                   expansion, scale, scale_exp
  USE Spectra, ONLY: DS1F, DS2F, DS3F, DS1H, DS2H, DS3H,   & 
                     exp_min, exp_range, l_val, freq_bins, &
                     Xo1, Xo2, Xo3, Yo1, Yo2, Yo3,         &
                     fber, erfb, loglem,                   &
                     time_window, time_bins, time_resolution
  USE Wind, ONLY: Average_kappa, SWS, nos, Ne_1AU, Ne_1Rs, Te_1AU, & 
                  Ti_1AU, Bo, Vsw_1AU
  USE Beam, ONLY: vel, VR, dv, halfdv, vmaxonVR 
  IMPLICIT NONE
  INTEGER :: i, j, io_error, alloc_error
  CHARACTER(LEN=8) :: DATE
  CHARACTER(LEN=10) :: TIME
  LOGICAL :: diditopen
!
333 FORMAT('START DATE: ',A2,'/',A2,'/',A4,'     ', & 
           'START TIME: ',A2,':',A2,':',A2,/)
313 FORMAT('***********************************',/,                     &
           '*** RESOLUTIONS/ARRAY SIZES ETC ***',/,                     &
           '***********************************',//,                    &
           '(DS)=Dynamic Spectra  ;  (GS)=Global Shock  ;  (R)=Ripple', & 
           //, 'Number of Frequency Bins (DS)=',I4,//,                  & 
           'Number of Time Bins (DS) in time window=',I4,//,            & 
           'Time Window (DS) kept in memory=',I5,' s', //,              & 
           'Calculation time resolution=',I5,' s', //,                  & 
           'Average Number of Ripples Packing Either Side of the Nose (GS)=',F5.1,//, & 
           'Shock Array Size (R)=',I6,//,                               & 
           'Calculation Region (R) Array Size x=',I4,//,                & 
           'Calculation Region (R) Array Size R=',I4,//,                & 
           'Electron Distribution Function Array Size V=',I6,//,        & 
           'Heliocentric Parameter Array Size P=',I4,///)   
370 FORMAT('*********************************',/,  &
           '*** PARAMETER RANGES of PLOTS ***',/,  &
           '*********************************',//, &
           'Maxtheta=', F13.8,/,                   &  
           'Mintheta=', F13.8,//,                  & 
           'MaxU=', ES10.1,/,                      & 
           'MinU=', ES10.1,//,                     & 
           'Maxb=', F6.1,/,                        & 
           'Minb=', F7.1,//,                       & 
           'MaxB1=', F6.1,/,                       & 
           'MinB1=', F6.1,//,                      & 
           'MaxNe=', F6.1,/,                       & 
           'MinNe=', F5.1,//,                      & 
           'MaxTe=', ES7.1,/,                      & 
           'MinTe=', ES7.1,//,                     & 
           'MaxTi=', ES7.1,/,                      &  
           'MinTi=', ES7.1,///)
!
  CALL DATE_AND_TIME(DATE,TIME)
!
  UT_time=0
  t_step=0
  first_pass=1
  DS1F=0.0_L ; DS1H=0.0_L
  DS2F=0.0_L ; DS2H=0.0_L
  DS3F=0.0_L ; DS3H=0.0_L
!   
  DO i=1,VR
    vel(i)=i*vmaxonVR 
  END DO    
! calculates parallel velocity spanned by one array point: 
  dv=vel(2)-vel(1)
  halfdv=0.5_L*dv
! 
! User input:
!
  OPEN(unit=20,file='./Dynamic_Spectra.input',status='old', & 
       IOSTAT=io_error)
    IF(io_error == 0)THEN 
      READ(20,'(17x,A)')Output_dir
      READ(20,'(11x,A)')file_ident
      READ(20,'(14x,F5.1)')Radial_Height
      READ(20,'(11x,F5.1)')Shock_Stop
      READ(20,'(5x,E8.1)')UCME
      READ(20,'(10x,E6.1)')CME_accel
      READ(20,'(10x,F4.1)')expansion
      READ(20,'(6x,F4.1)')scale
      READ(20,'(10x,F4.2)')scale_exp
      READ(20,'(8x,E7.1)')exp_min
      READ(20,'(10x,E6.1)')exp_range
      READ(20,'(6x,F3.1)')l_val
      READ(20,'(14x,F3.1)')Average_kappa
      READ(20,'(7x,E5.1)')Ne_1AU
      READ(20,'(7x,E6.2)')Ne_1Rs
      READ(20,'(7x,E5.1)')Te_1AU
      READ(20,'(7x,E5.1)')Ti_1AU
      READ(20,'(3x,E13.8)')Bo
      READ(20,'(8x,E5.1)')Vsw_1AU
      READ(20,'(4x,E10.5)')Xo1
      READ(20,'(4x,E10.5)')Yo1
      READ(20,'(4x,E10.5)')Xo2
      READ(20,'(4x,E10.5)')Yo2
      READ(20,'(4x,E10.5)')Xo3
      READ(20,'(4x,E10.5)')Yo3
      READ(20,'(9x,F13.8)')Maxtheta
      READ(20,'(9x,F13.8)')Mintheta
      READ(20,'(5x,E10.1)')MaxU
      READ(20,'(5x,E10.1)')MinU
      READ(20,'(5x,F6.1)')Maxb
      READ(20,'(5x,F7.1)')Minb
      READ(20,'(6x,F6.1)')MaxB1
      READ(20,'(6x,F6.1)')MinB1
      READ(20,'(6x,F6.1)')MaxNe
      READ(20,'(6x,F5.1)')MinNe
      READ(20,'(6x,E7.1)')MaxTe
      READ(20,'(6x,E7.1)')MinTe
      READ(20,'(6x,E7.1)')MaxTi
      READ(20,'(6x,E7.1)')MinTi    
      READ(20,'(4x,I1)')nos
      ALLOCATE(SWS(nos), STAT=alloc_error)
      IF(alloc_error /= 0)THEN 
        print*,'Failed to ALLOCATE array of SWS(',nos,')'
        STOP
      END IF
      DO i= 1,nos
        READ(20,'(13x,I1)')SWS(i)%shape
        READ(20,'(11x,I1)')SWS(i)%nop
        ALLOCATE(SWS(i)%param(SWS(i)%nop), STAT=alloc_error)
        IF(alloc_error /= 0)THEN 
          print*,'Failed to ALLOCATE array of SWS(',i, & 
                 ')%param(',SWS(i)%nop,')'
          STOP
        END IF
        DO j=1,SWS(i)%nop
          READ(20,'(16x,E8.2)')SWS(i)%param(j)
        END DO
        READ(20,'(6x,E8.2)')SWS(i)%theta 
        READ(20,'(3x,E8.2)')SWS(i)%Ne
        READ(20,'(3x,E8.2)')SWS(i)%B1
        READ(20,'(3x,E8.2)')SWS(i)%Ti
        READ(20,'(3x,E8.2)')SWS(i)%Te
        READ(20,'(4x,E8.2)')SWS(i)%Vsw
        READ(20,'(2x,E8.2)')SWS(i)%b
        READ(20,'(6x,E8.2)')SWS(i)%Kappa
      END DO
    ELSE 
      print*, 'No Dynamic_Spectra.input file found,' 
      print*, 'proceeding using default values...' 
      Output_dir='data/spectra/'
      file_ident='temp.dat' 
      Radial_Height=1.1 
      Shock_Stop=214.78679  !1AU
      UCME=1500.0e3 
      CME_accel=5.9
      expansion=0.0
      scale=1.0
      scale_exp=0.37
      exp_min=4.0
      exp_range=5.0
      l_val=1.0
      Average_kappa=3.0
      Ne_1AU=5.1e6
      Ne_1Rs=1.68e14
      Te_1AU=2.4e5
      Ti_1AU=1.7e5
      Bo=1.89743128e-4
      Vsw_1AU=4.0e5
      Xo1=1.49599e11
      Yo1=1.00e9
      Xo2=1.00e11
      Yo2=1.00e11
      Xo3=1.00e11
      Yo3=-1.00e11
      Maxtheta=6.283185307
      Mintheta=3.141592654
      MaxU=150.0e4
      MinU=0.0e1
      Maxb=-7.0
      Minb=-10.0
      MaxB1=-3.0
      MinB1=-9.0
      MaxNe=15.0
      MinNe=6.0
      MaxTe=2.0e6 
      MinTe=4.0e4
      MaxTi=3.0e6 
      MinTi=4.0e4
      nos=0
    END IF 
  INQUIRE(unit=20,EXIST=diditopen)
  IF(diditopen)CLOSE(unit=20)
!
  Radial_Height=Radial_Height*Sol_rad
  Shock_Stop=Shock_Stop*Sol_rad
!
  fber=freq_bins/exp_range
  erfb=exp_range/freq_bins
  loglem=LOG10(l_val)+exp_min
!
  OPEN(unit=21,file=TRIM(Output_dir)//'LOG_'//file_ident,status='new')
     write(unit=21,fmt=333)DATE(7:8),DATE(5:6),DATE(1:4), & 
                           TIME(1:2),TIME(3:4),TIME(5:6)
     write(unit=21,fmt=*)'HUGE(REAL(KIND=L))=',HUGE(MaxU)
     write(unit=21,fmt=313)freq_bins,time_bins,time_window,      &
                           time_resolution,Avg_num_per_timestep, & 
                           SR, xR, RR, VR, PR
     write(unit=21,fmt=370)Maxtheta, Mintheta, MaxU, MinU, Maxb, Minb, & 
                           MaxB1, MinB1, MaxNe, MinNe, MaxTe, MinTe,   & 
                           MaxTi, MinTi
  ENDFILE(unit=21)
  CLOSE(unit=21)           
!
  END SUBROUTINE Initialise_Global
!____________________________________________________________________

  SUBROUTINE Initialise_Time_Step
! *******************************
  USE Global, ONLY: L, AU, Sol_rad, Sol_period, pi
  USE Shock, ONLY: rip_pck_num, Radial_Height, SX, SY, UT_time, t_step, & 
                   UCME, CME_accel, CME_b, scale, scale_exp
  USE Wind, ONLY:  Average_Vsw, Vsw_1AU, Bo, Average_B1
  USE Spectra, ONLY: spec_shift, time_window, time_bins, time_resolution
  USE OutPut, ONLY: Output_dir, file_ident,    & 
                    theta_arr, U_arr, b_arr, B1_arr, & 
                    Ne_arr, Te_arr, Ti_arr
  IMPLICIT NONE
!
  rip_pck_num=0.5_L ! 0.5 cause calcing 2 rips per circ on global
  UT_time=UT_time+t_step
!
  UCME=UCME-CME_accel*t_step 
!
  Radial_Height=Radial_Height+(UCME*t_step)
  SX=0.1e-7_L ; SY=0.1e-7_L  ! SY cannot be zero (see m=)
                             ! could probably solve this by requiring 
                             ! observer to never be on Y=0 line
! 
! Calculate solar wind and shock parameters with heliocentric variations:
!
  CME_b=(scale/Radial_Height)/(Radial_Height/AU)**scale_exp 
  Average_Vsw=Vsw_1AU*(Radial_Height/AU)**0.19  ! R&C 1998 III
  Average_B1=SQRT((Bo*Sol_rad*Sol_rad/(Radial_Height*Radial_Height))**2 & 
                    +(Bo*(2.0*pi/Sol_period)*Sol_rad*Sol_rad            & 
                    /(Average_Vsw*Radial_Height))**2)
!
! calculate time step between shock calculations based on correlation
! time; correlation time approximated by time for shock to traverse 
! a distance corresponding to the scale size of a ripple:
!
  t_step=1.37e5/(SQRT(Average_B1)*(UCME-Average_Vsw))
  IF(t_step <= time_resolution)t_step=time_resolution 
!
  IF(t_step >= time_window)THEN 
    OPEN(unit=21, file=TRIM(Output_dir)//'LOG_'//file_ident, & 
         status='old', position='append')
      write(unit=21,fmt=*) 'Whoops!!!' 
      write(unit=21,fmt=*) 't_step=',t_step 
      write(unit=21,fmt=*) 'Radial_Height=',Radial_Height 
    ENDFILE(unit=21) 
    CLOSE(unit=21) 
  END IF
! 
! How much(in time) of the dynamic spectra window to dump to a file: 
!
  spec_shift=t_step*time_bins/(1.0_L*time_window)
! 
  t_step=spec_shift*((1.0_L*time_window)/time_bins)
!
! Initialise Range of parameter fluctuations: 
!  
  theta_arr=0 ; U_arr=0 ; b_arr=0 ; B1_arr=0
  Ne_arr=0 ; Te_arr=0 ; Ti_arr=0
! 
  END SUBROUTINE Initialise_Time_Step
!____________________________________________________________________

  SUBROUTINE Ripple_Parameters (WHICH_SIDE)
! *****************************************  
  USE Global, ONLY: L, c, csq, kbolt, AU, Sol_rad, Sol_period, mi, me, & 
                    pi, rtpi, sqthree, calculate
  USE Wind, ONLY: Average_theta, Average_B1, Average_Ti, Average_Ne,  & 
                  Average_Te, Average_Vsw, Te, Ti, Ne, B1, Vsw, Vd,   & 
                  theta, Ve, Vesq, cs, kappa, Average_kappa, stheta,  & 
                  stheta2, ctheta, tantheta, stheta1on, theta_Parker, & 
                  nos, Average_fp, Te_1AU, Ti_1AU, Bo, Ne_1AU, Ne_1Rs, Vsw_1AU
  USE OutPut, ONLY: theta_arr, U_arr, b_arr, B1_arr, Ne_arr, Te_arr,  & 
                    Ti_arr, PR, Maxtheta, Mintheta, MaxU, MinU, Maxb, & 
                    Minb, MaxB1, MinB1, MaxNe, MinNe, MaxTe, MinTe,   & 
                    MaxTi, MinTi
  USE Emis, ONLY: temp4gammaLS, temp4zetaF, erftemp1, erftemp2, & 
                  temp4phiF, temp4phiH, Vcrit, gmemisq
  USE Shock, ONLY: Radial_Height, Ripple_Height, RHonAU, b, Average_b, & 
                   CME_b, U, UCME, SX, SY, SY_temp, theta_radial,      & 
                   alpha_global, expansion, RHSX
  USE Beam, ONLY: b4on, bsqsthetasq, bstheta1on, Uctp
  USE NR, ONLY: Gaussian_Random
  IMPLICIT NONE 
  INTEGER, INTENT(IN) :: WHICH_SIDE            
  REAL(KIND=L),PARAMETER:: Uc=1.8_L            
  REAL(KIND=L) :: gamma, migmesq, Vecu
  REAL(KIND=L), DIMENSION(10) :: BOB=0.0_L 
  INTEGER :: i
!
  RHSX=Radial_Height-SX  !X-position of ripple wrt sun 
  SY_temp=((-1)**WHICH_SIDE)*SY   !Y-position of ripple wrt sun
!
  CALL Line_of_Sight
!
  Ripple_Height=SQRT(RHSX*RHSX + SY*SY)
  RHonAU=Ripple_Height/AU
!
  Average_Te=Te_1AU/RHonAU**(0.42) ! R&C 1998 III
  Average_Ti=Ti_1AU/RHonAU**(0.67)  !nominal 1AU gives coronal ~1.8 MK 
  Average_Ne=(Ne_1AU/RHonAU**2.19                     &  ! R&C 1998 III
             +Ne_1Rs*(Ripple_Height/Sol_rad)**(-6.13))   ! SAITO POLAND MUNRO 1977
!
!Calculate Vsw again for each ripple (specific radial height)
  Average_Vsw=Vsw_1AU*(RHonAU)**0.19  ! R&C 1998 III
!
   Average_fp=8.976_L*SQRT(Average_Ne)    
!  
! Theta and B1 calculations based on Parker Spiral:
!
  theta_Parker=ATAN(2.0*pi*Ripple_Height/(Sol_period*Average_Vsw) & 
                    *(1.0-Sol_rad/Ripple_Height)) 
  theta_radial=ATAN(SY/RHSX)
  alpha_global=ATAN(2.0*CME_b*SY)
!
  IF(WHICH_SIDE == 1)Average_theta=theta_Parker-alpha_global+theta_radial
  IF(WHICH_SIDE == 2)Average_theta=theta_Parker+alpha_global-theta_radial
  IF(Average_theta < 0.0)Average_theta=-Average_theta
!
  Average_B1=SQRT((Bo*Sol_rad*Sol_rad/(Ripple_Height*Ripple_Height))**2 & 
                    +(Bo*(2.0*pi/Sol_period)*Sol_rad*Sol_rad            & 
                    /(Average_Vsw*Ripple_Height))**2)
!
  Average_b=SQRT(Average_B1)*7.29927e-6 !Mancuso&Spangler 1999 
                                        !(ApJ 525, pg195)
                                        !and references therein
!
!  BOB(1)=+-50% ; BOB(2:4)=+-10% ; BOB(5:6)=+-20% ; BOB(7:10)=+-30%
!
  CALL Gaussian_Random(BOB)  
!
  WHERE(BOB < -1.0)
    BOB=-1.0_L 
  END WHERE
  WHERE(BOB > 1.0) 
    BOB=1.0_L
  END WHERE
!
  BOB(1)=1.0_L+BOB(1)*0.5_L      
  BOB(2:4)=1.0_L+BOB(2:4)*0.1_L  
  BOB(5:6)=1.0_L+BOB(5:6)*0.2_L  
  BOB(7:10)=1.0_L+BOB(7:10)*0.3_L
!  
  theta=2.0*pi-BOB(1)*Average_theta 
  Ne=BOB(6)*Average_Ne
  B1=BOB(4)*Average_B1
  Ti=BOB(10)*Average_Ti
  Te=BOB(8)*Average_Te
  Vsw=BOB(9)*Average_Vsw 
  b=BOB(5)*Average_b
  Kappa=BOB(7)*Average_Kappa 
!
  IF(nos /= 0)CALL Structures
!
  stheta=SIN(theta)
  stheta2=2.0_L*stheta
  stheta1on=1.0_L/stheta
  ctheta=COS(theta)
  tantheta=TAN(theta)
  U=UCME*(COS(alpha_global)+(expansion-expansion*COS(alpha_global))) & 
    -Vsw*COS(theta_radial- alpha_global)  
!    
  Vd=-U*stheta
  Ve=SQRT((kbolt/me)*Te)
  cs=SQRT(5.0_L/3.0_L*kbolt*(Te+Ti)/mi)
!
! bin selection is
! ArraySize * (Value-MinOfRange)/Range
! 
  i=PR*(theta-Mintheta)/(Maxtheta-Mintheta)
  IF(i < 1)i=1 ; IF(i > PR)i=PR
  theta_arr(i)=theta_arr(i)+1
!
  i=PR*(U-MinU)/(MaxU-MinU)
  IF(i < 1)i=1 ; IF(i > PR)i=PR
  U_arr(i)=U_arr(i)+1 
!
  i=PR*(LOG10(b)-Minb)/(Maxb-Minb)
  IF(i < 1)i=1 ; IF(i > PR)i=PR
  b_arr(i)=b_arr(i)+1 
!
  i=PR*(LOG10(B1)-MinB1)/(MaxB1-MinB1)
  IF(i < 1)i=1 ; IF(i > PR)i=PR
  B1_arr(i)=B1_arr(i)+1 
!
  i=PR*(LOG10(Ne)-MinNe)/(MaxNe-MinNe)
  IF(i < 1)i=1 ; IF(i > PR)i=PR
  Ne_arr(i)=Ne_arr(i)+1
!
  i=PR*(Te-MinTe)/(MaxTe-MinTe) 
  IF(i < 1)i=1 ; IF(i > PR)i=PR
  Te_arr(i)=Te_arr(i)+1 
!
  i=PR*(Ti-MinTi)/(MaxTi-MinTi)
  IF(i < 1)i=1 ; IF(i > PR)i=PR
  Ti_arr(i)=Ti_arr(i)+1
! 
  IF(theta > 5.8)THEN !< 30deg of shock normal or > 60deg from perp
    calculate=0       !was (Average_theta < 0.5) until 23/6/2003
    RETURN 
  END IF
!
  Vesq=Ve*Ve
  Vecu=Ve*Ve*Ve
  gamma=1.0_L+3.0_L*Ti/Te 
  migmesq=SQRT(mi/(gamma*me))
  gmemisq=SQRT((gamma*me)/mi)
  Vcrit=1.5_L*migmesq*Ve
  temp4gammaLS=4.0_L*migmesq*Vesq
  temp4zetaF=(4.0_L*gamma*me)/(45.0_L*mi)
  erftemp1=(Ve*sqthree/c)+(2.0_L/3.0_L)*gmemisq
  erftemp2=(Ve*sqthree/c)-(2.0_L/3.0_L)*gmemisq 
  temp4phiF=72.0_L*EXP(-Uc*Uc)/(Uc*rtpi)*sqthree*Vecu/(csq*c)
  temp4phiH=18.0_L*sqthree*migmesq/(5.0_L*gamma*csq*c)*(Vecu/csq)
!
  b4on=4.0_L/b
  bsqsthetasq=1.0_L/(b*b *stheta*stheta)
  bstheta1on=1.0_L/(b *stheta)
  Uctp=U*COS(theta-pi)
!
  CALL bckgrnd_distribution
!
  CALL shock_parameter_arrays
!
  END SUBROUTINE Ripple_Parameters 
!____________________________________________________________________

  SUBROUTINE Line_of_Sight
! ************************
  USE Global, ONLY: L
  USE Shock, ONLY: CME_b, SY_temp, Radial_Height, RHSX
  USE Spectra, ONLY: Xo1,Xo2,Xo3,Yo1,Yo2,Yo3,blocked1,blocked2,blocked3
  IMPLICIT NONE
  REAL(KIND=L) :: m, Y1, Y2
  REAL(KIND=L) :: CME_b2, CME_b4
!
  CME_b2=2.0_L*CME_b
  CME_b4=4.0_L*CME_b
!  
!  determine whether line of sight is clear or passes through shock and 
!  will thus need frequency filtering of final DS
!
  m= (Xo1-RHSX)/(Yo1-SY_temp)
  Y1=(SQRT(m*m+CME_b4*(Radial_Height+m*Yo1-Xo1))-m)/CME_b2
  Y2=-(SQRT(m*m+CME_b4*(Radial_Height+m*Yo1-Xo1))+m)/CME_b2
!
  IF(ABS(SY_temp-Y1) > ABS(SY_temp-Y2))THEN !soltn Y1orY2 isn't SY_temp
    IF(ABS(Yo1-SY_temp) > ABS(Yo1-Y1))THEN  !is second soln on the shock  
      blocked1= 1                           !closer to 
    ELSE                                    !or 
      blocked1= 0                           !further from observer 
    END IF 
  ELSE
    IF(ABS(Yo1-SY_temp) > ABS(Yo1-Y2))THEN !is second soln on the shock 
      blocked1= 1                          !closer to 
    ELSE                                   !or 
      blocked1= 0                          !further from observer 
    END IF 
  END IF
!
  m= (Xo2-RHSX)/(Yo2-SY_temp)
  Y1=(SQRT(m*m+CME_b4*(Radial_Height+m*Yo2-Xo2))-m)/CME_b2
  Y2=-(SQRT(m*m+CME_b4*(Radial_Height+m*Yo2-Xo2))+m)/CME_b2
!
  IF(ABS(SY_temp-Y1) > ABS(SY_temp-Y2))THEN !soltn Y1orY2 isn't SY_temp
    IF(ABS(Yo2-SY_temp) > ABS(Yo2-Y1))THEN  !is second soln on the shock 
      blocked2= 1                           !closer to 
    ELSE                                    !or 
      blocked2= 0                           !further from observer 
    END IF 
  ELSE
    IF(ABS(Yo2-SY_temp) > ABS(Yo2-Y2))THEN !is second soln on the shock 
      blocked2= 1                          !closer to  
    ELSE                                   !or 
      blocked2= 0                          !further from observer 
    END IF 
  END IF  
!
  m= (Xo3-RHSX)/(Yo3-SY_temp)
  Y1=(SQRT(m*m+CME_b4*(Radial_Height+m*Yo3-Xo3))-m)/CME_b2
  Y2=-(SQRT(m*m+CME_b4*(Radial_Height+m*Yo3-Xo3))+m)/CME_b2
!
  IF(ABS(SY_temp-Y1) > ABS(SY_temp-Y2))THEN !soltn Y1orY2 isn't SY_temp
    IF(ABS(Yo3-SY_temp) > ABS(Yo3-Y1))THEN  !is second soln on the shock 
      blocked3= 1                           !closer to 
    ELSE                                    !or 
      blocked3= 0                           !further from observer 
    END IF 
  ELSE
    IF(ABS(Yo3-SY_temp) > ABS(Yo3-Y2))THEN !is second soln on the shock 
      blocked3= 1                          !closer to 
    ELSE                                   !or 
      blocked3= 0                          !further from observer 
    END IF 
  END IF
! 
  END SUBROUTINE Line_of_Sight
!____________________________________________________________________

  SUBROUTINE Structures
! *********************
  USE Global, ONLY: L, pi, Sol_rad
  USE Wind, ONLY: SWS, nos, theta, Ne, B1, Te, Ti, Vsw, Kappa
  USE Shock, ONLY: SY_temp, Ripple_Height, b, RHSX
  IMPLICIT NONE
  REAL(KIND=L) :: R_spiral
  INTEGER :: j
!  
  DO j=1,nos
    SELECT CASE (SWS(j)%shape)
      CASE (1) !CIR  
        R_spiral=137*(SWS(j)%param(1)*pi-ATAN(SY_temp/RHSX))*Sol_rad
        IF(Ripple_Height > 7*Sol_rad)THEN
          IF(((SY_temp >= 0) .AND. (RHSX <= ((137*SWS(j)%param(1)*pi*  & 
            Sol_rad)*(1+SWS(j)%param(2))))) .OR. ((SY_temp <  0) .AND. &
            (RHSX > ((137*SWS(j)%param(1)*pi*Sol_rad)*                 & 
            (1-2*SWS(j)%param(2))))))THEN
            IF(((R_spiral-Ripple_Height) <= 0) .AND. ((R_spiral*   & 
               (1+SWS(j)%param(2))-Ripple_Height) > 0))THEN  !!!HIGH 
              theta=theta*SWS(j)%theta
              Ne=Ne*(1-SWS(j)%Ne)
              B1=B1*(1-SWS(j)%B1)
              Ti=Ti*(1-SWS(j)%Ti)
              Te=Te*(1+SWS(j)%Te)
              Vsw=Vsw*(1+SWS(j)%Vsw)
              b=b*SWS(j)%b
              Kappa=Kappa*(1+SWS(j)%Kappa) 
            END IF  
            IF(((R_spiral-Ripple_Height) > 0) .AND. ((R_spiral* & 
                (1-SWS(j)%param(2))-Ripple_Height) <= 0))THEN  !!!MID
              theta=theta*SWS(j)%theta
              Ne=Ne*(1+SWS(j)%Ne)
              B1=B1   !*(1+SWS(j)%B1)
              Ti=Ti*(1+3*SWS(j)%Ti)
              Te=Te*(1+3*SWS(j)%Te)
              Vsw=Vsw*(1+2*SWS(j)%Vsw)
              b=b*SWS(j)%b
              Kappa=Kappa*(1-SWS(j)%Kappa)  
            END IF  
            IF(((R_spiral*(1-SWS(j)%param(2))-Ripple_Height) > 0) .AND. & 
              ((R_spiral*(1-2*SWS(j)%param(2))-Ripple_Height) < 0))THEN  !!!LOW
              theta=theta*SWS(j)%theta
              Ne=Ne*(1+SWS(j)%Ne)
              B1=B1*(1+6*SWS(j)%B1)
              Ti=Ti*(1+10*SWS(j)%Ti)
              Te=Te*(1+10*SWS(j)%Te)
              Vsw=Vsw*(1+8*SWS(j)%Vsw)
              b=b*SWS(j)%b
              Kappa=Kappa*(1-SWS(j)%Kappa) 
            END IF  
          END IF 
        END IF
      CASE (2)  !annulus between heliocentric radius aaa & bbb
        IF((((RHSX*RHSX)+(SY_temp*SY_temp)) > (SWS(j)%param(1)*    & 
           SWS(j)%param(1))).AND. (((RHSX*RHSX)+(SY_temp*SY_temp)) & 
            < (SWS(j)%param(2)*SWS(j)%param(2))))THEN 
          theta=theta*SWS(j)%theta
          Ne=Ne*SWS(j)%Ne
          B1=B1*SWS(j)%B1
          Ti=Ti*SWS(j)%Ti
          Te=Te*SWS(j)%Te
          Vsw=Vsw*SWS(j)%Vsw
          b=b*SWS(j)%b
          Kappa=Kappa*SWS(j)%Kappa  
        END IF    
      CASE (4)   !square/rectangular region bound by aaa, bbb, ccc, ddd
        IF((RHSX > SWS(j)%param(1)) .AND. (RHSX < SWS(j)%param(2)) &
           .AND. (SY_temp < SWS(j)%param(3)) .AND.                 &    
           (SY_temp > SWS(j)%param(4)))THEN  
          theta=theta*SWS(j)%theta
          Ne=Ne*SWS(j)%Ne
          B1=B1*SWS(j)%B1
          Ti=Ti*SWS(j)%Ti
          Te=Te*SWS(j)%Te
          Vsw=Vsw*SWS(j)%Vsw
          b=b*SWS(j)%b
          Kappa=Kappa*SWS(j)%Kappa   
        END IF
      CASE (3) !circular region of radius aaa, & center (X,Y)=(bbb,ccc)
        IF(((RHSX-SWS(j)%param(2))**2+(SY_temp-SWS(j)%param(3))**2) & 
            < (SWS(j)%param(1)*SWS(j)%param(1)))THEN  
          theta=theta*SWS(j)%theta
          Ne=Ne*SWS(j)%Ne
          B1=B1*SWS(j)%B1
          Ti=Ti*SWS(j)%Ti
          Te=Te*SWS(j)%Te
          Vsw=Vsw*SWS(j)%Vsw
          b=b*SWS(j)%b
          Kappa=Kappa*SWS(j)%Kappa   
        END IF
      CASE (6) !coronal loop
        IF((((RHSX-SWS(j)%param(2))**2+(SY_temp-SWS(j)%param(3))**2) & 
           > (SWS(j)%param(1)*SWS(j)%param(1)))  .AND.               & 
           (((RHSX-SWS(j)%param(5))**2+(SY_temp-SWS(j)%param(6))**2) & 
           < (SWS(j)%param(4)*SWS(j)%param(4))))THEN  
          IF(Ripple_Height > 0.75*SWS(j)%param(1)+Sol_rad)THEN 
            theta=theta+0.25*pi   !*SWS(j)%theta !4.95 
          ELSE IF(Ripple_Height > 0.5*SWS(j)%param(1)+Sol_rad)THEN  
            theta=theta+0.125*pi 
          END IF 
          Ne=Ne*SWS(j)%Ne    !5.0e15
          B1=B1*SWS(j)%B1    !0.75e-7 
          Ti=Ti*SWS(j)%Ti    !7.5e4
          Te=Te*SWS(j)%Te    !7.5e4
          Vsw=SWS(j)%Vsw     !0.0
          b=b*SWS(j)%b       !b
          Kappa=Kappa*SWS(j)%Kappa  
        END IF
    END SELECT
  END DO
!
  END SUBROUTINE Structures
!____________________________________________________________________

  SUBROUTINE bckgrnd_distribution 
! *******************************    
  USE Global, ONLY: L, csq, rtpi
  USE Beam, ONLY: gammafn, ConstDist, vel, vparray, temp4vponvd
  USE Wind, ONLY: Kappa, Ve, Vd, Vesq
  IMPLICIT NONE
!
! For parallel velocity from 0 to vmax, calculate reduced
! solar wind distribution function.
!
  ConstDist=(gammafn(Kappa+1.0_L)/(Kappa*rtpi*gammafn(Kappa-0.5_L))  & 
            *Ve**(2.0_L*Kappa-1.0_L))*1.0e-6**Kappa
! 
  vparray=ConstDist*(((Vesq+vel*vel)*1.0e-6)**(-Kappa) & 
          -((Vesq+vel*vel+csq)*1.0e-6)**(-Kappa))
!  
  temp4vponvd=vel/Vd
!  
  END SUBROUTINE bckgrnd_distribution
!____________________________________________________________________
 
  SUBROUTINE shock_parameter_arrays
! *********************************  
  USE Global, ONLY: L, pi, pi15, mu0, mi, calculate
  USE Shock, ONLY: b, Rrange, xrange, Srange, SRonSrange, SR, B2B1,  & 
                   UnSalphonB2B1, alph, thetabn, Vc, Un, U, UnSalph, & 
                   Ma, Mssq, S_limit
  USE Wind, ONLY: stheta, ctheta, tantheta, stheta1on, Vd,    & 
                  B1, Ne, cs, theta
  IMPLICIT NONE
!
  REAL(KIND=L) :: Rs, xs, bsthetactheta2, ctheta1on, s2theta
  INTEGER :: S
  REAL(KIND=L) :: Va   ! Alfven speed
!
  Va=B1/SQRT(mu0*mi*Ne)
!
  bsthetactheta2=2.0_L*b*stheta*ctheta
  ctheta1on=1.0_L/ctheta
  s2theta=SIN(2.0_L*theta)
  Rrange=3.0_L*(64/(LOG10(b))**2)/b  !64=8^2   approx emission  coverage
  xrange=0.3_L*(100/(LOG10(b))**2)/b !100=10^2 approx emission  coverage
!  
  Srange=(b*xrange+(1.0_L/tantheta) &
        *SQRT(-b*xrange*tantheta/ctheta**3))/(b*tantheta)
  SRonSrange=SR/Srange
!
 S_limit=0
  DO S=1,SR
    Rs=S*Srange/SR
    xs=-ctheta1on*ctheta1on*(stheta1on+SQRT(-4.0_L*b*Rs*ctheta &
       +stheta1on*stheta1on)-b*Rs*s2theta)/(2.0_L*b)
!
!   For position RS,xs on shock calculate the angle(alph)
!   between B-field and shock tangent and the cutoff
!   velocity(Vc): 
!      
    alph(S)=ATAN(tantheta/(1.0_L+1.0_L/(bsthetactheta2  &
            *(-Rs*stheta+xs*ctheta))))
    thetabn(S)=pi*0.5_L-alph(S)
    Vc(S)= ABS(Vd*(1.0_L/TAN(alph(S))))
    Un(S)= U*COS(pi15-theta+alph(S))  ! angle between tangent to shock
!                                   surface and tangential B-field line;
!
    UnSalph(S)=Un(S)/SIN(alph(S))
!
    IF(alph(S) > 1.047)THEN !~60 degrees (59.988681...)
      S_limit=S
      EXIT
    END IF
!
  END DO
  IF(S_limit==0)S_limit=SR
!
  Ma(1:S_limit)=Un(1:S_limit)/Va
  Mssq(1:S_limit)=(Un(1:S_limit)*Un(1:S_limit))/(cs*cs)
!
  IF(MAXVAL(Ma(1:S_limit)) <= 1.2)THEN !Don't bother calculating weak shock
    calculate=0
    RETURN 
  END IF
!
  CALL B2B1_ephi
  UnSalphonB2B1(1:S_limit)=UnSalph(1:S_limit)*UnSalph(1:S_limit) & 
                           /B2B1(1:S_limit)
!
  END SUBROUTINE shock_parameter_arrays
!____________________________________________________________________

  SUBROUTINE B2B1_ephi
! *********************                                               
  USE Global, ONLY: L, kbolt, me2on
  USE Shock, ONLY: Ma, n2n1, B2B1, B2B1_1, thetabn, ephi, ephi2onme, & 
                   S_limit
  USE Wind, ONLY: Te
  IMPLICIT NONE
!
  REAL(KIND=L), ALLOCATABLE, DIMENSION(:) :: Bt2Bt1    ! tang comp B2/B1
  REAL(KIND=L), ALLOCATABLE, DIMENSION(:) :: cthetabn  ! Cos(thetabn)
  REAL(KIND=L), ALLOCATABLE, DIMENSION(:) :: sthetabn  ! Sin(thetabn)
  REAL(KIND=L), ALLOCATABLE, DIMENSION(:) :: ephipara
  REAL(KIND=L), ALLOCATABLE, DIMENSION(:) :: ephiperp
  REAL(KIND=L), PARAMETER:: A=1.1_L   ! temperature anisotropy
  REAL(KIND=L), PARAMETER :: zeta=-0.1_L  ! anisotropy jump
  REAL(KIND=L) :: kTeper, kTepar
  INTEGER :: Alloc_error, Dealloc_error
!
  CALL cgroots
!
  ALLOCATE(Bt2Bt1(S_limit), cthetabn(S_limit), sthetabn(S_limit), &
           ephipara(S_limit), ephiperp(S_limit), STAT=alloc_error)
  IF(alloc_error /= 0)THEN 
    print*,'Failed to ALLOCATE arrays for Bt2Bt1, cthetabn, sthetabn, cotthetabn'
    STOP
  END IF
!
  cthetabn=COS(thetabn(1:S_limit))
  sthetabn=SIN(thetabn(1:S_limit))
  Bt2Bt1=n2n1(1:S_limit)*(Ma(1:S_limit)*Ma(1:S_limit)-cthetabn*cthetabn) & 
         /(Ma(1:S_limit)*Ma(1:S_limit)-n2n1(1:S_limit)*cthetabn*cthetabn)
  B2B1(1:S_limit)=SQRT(cthetabn*cthetabn+sthetabn*sthetabn*Bt2Bt1*Bt2Bt1)
  B2B1_1(1:S_limit)=B2B1(1:S_limit)-1.0_L
  kTeper=(3.0_L*Te /(A+2.0_L))*kbolt ;  kTepar=A*kTeper
!
  ephi(1:S_limit)=0.0_L
  WHERE(B2B1(1:S_limit) > 1.0)
    ephiperp=kTeper*(B2B1(1:S_limit)-1.0_L)
!   The following model assumes n(x) propto B(x):
    ephipara=kTepar*(1.0_L + zeta)*(B2B1(1:S_limit)-1.0_L)
    ephi(1:S_limit)=ephipara+ephiperp
  END WHERE
  WHERE(ephi(1:S_limit) == 0.0_L)ephi(1:S_limit)=MAXVAL(ephi(1:S_limit))
!
  ephi2onme(1:S_limit)=me2on*ephi(1:S_limit)
!
  DEALLOCATE(Bt2Bt1, cthetabn, sthetabn, ephipara, ephiperp, & 
             STAT=dealloc_error)
  IF(dealloc_error /= 0)PRINT*,'Deallocation error for Bt2Bt1 & thetabn arrays'
!
  END SUBROUTINE B2B1_ephi
!____________________________________________________________________

  SUBROUTINE cgroots 
! *******************
!
!   Uses a Numerical Recipes routine, ZROOTS, to calculate
!   all three complex roots of the cubic equation for the 
!   density jump across a shock as a function of theta,
!   M_a and M_s. The real part of one of the solutions
!   gives the density ratio n2/n1.
!
!   START DATE 5 May 1995 (Iver Cairns)
!   MODIFIED 11 Sept 2000 (Zdenka Kuncic) to incorporate
!                         into Fortran 90 code Foreshock
!                                               
  USE Global, ONLY:  L, pi
  USE Shock, ONLY: S_limit, Ma, Mssq, thetabn, n2n1
  IMPLICIT NONE  
  INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: SP = KIND(1.0) 
  INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
  INTEGER, PARAMETER :: LGT = KIND(.true.)
  INTEGER(I4B), PARAMETER :: M=3      ! Degree of the Polynomial
  COMPLEX(KIND=SPC) :: C_COEFF(M+1),C_ROOTS(M)
  REAL(KIND=L) :: gamma,GBETA
  REAL(KIND=L) :: GP1,GP2,GB,GM1,GM2
  REAL(KIND=L) ::REAL_ROOTS(M)
  REAL(KIND=L) :: COST,COS2,COS4,Ma6,Ma4,Ma2
  REAL(KIND=SP) :: A,BB,CC,D
  INTEGER :: S 
  LOGICAL(LGT) :: POLISH
!
    POLISH = .TRUE.      
!
    gamma = 5.0_L/3.0_L
    GP1 = gamma + 1._L ; GP2 = gamma + 2._L
    GM1 = gamma - 1._L ; GM2 = gamma - 2._L
!       
  DO S=1,S_limit
!
    IF((Ma(S) > 1.2) .AND. (Mssq(S) > 1.44))THEN
!
      COST = COS(pi-thetabn(S)) ; COS2 = COST*COST ; COS4 = COS2*COS2
      Ma2=Ma(S)*Ma(S) ; Ma4=Ma2*Ma2 ; Ma6=Ma4*Ma2
  !    
  ! X EQUATION: AX^3 + BX^2 + CX + D:
  !
      GBETA = 2.0_L * Ma2 / Mssq(S)
      GB = gamma + GBETA
      A = GP1 * Ma6
      BB = GM1*Ma6 + GP2*COS2*Ma4 + GB*Ma4 ; BB = -1.0_L * BB
      CC = (GM2 + gamma*COS2)*Ma4 + (2.*GBETA + GP1)*COS2*Ma2
      D = -1.0_L * (GM1*COS2*Ma2 + COS4 * GBETA)
  !
  ! Prepare the arrays for the root find:
  !
      C_COEFF(1) = CMPLX(D,0.0_sp)
      C_COEFF(2) = CMPLX(CC,0._sp)
      C_COEFF(3) = CMPLX(BB,0._sp)
      C_COEFF(4) = CMPLX(A,0._sp)
      C_ROOTS = (0.0_sp,0.0_sp)
  !               
  ! now find the roots:
  !
      CALL ZROOTS(C_COEFF,M,C_ROOTS,POLISH)
  !
      REAL_ROOTS(3)=REAL(C_ROOTS(3))  
  !   
      n2n1(S) = 1.0_L/REAL_ROOTS(3)
  !
    ELSE
      n2n1(S) = 1.0_L
    ENDIF
!
    IF(n2n1(S) < 1.0)n2n1(S) = 1.0_L
!
  END DO     
!
  END SUBROUTINE cgroots 
!____________________________________________________________________

 SUBROUTINE zroots(a,M,roots,polish)
!**********************************   
 IMPLICIT NONE 
 INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9) 
 INTEGER, PARAMETER :: SP = KIND(1.0) 
 INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
 INTEGER, PARAMETER :: LGT = KIND(.true.)
!
 INTEGER(I4B), INTENT(IN) :: M
 COMPLEX(SPC), DIMENSION(M+1), INTENT(IN) :: a 
 COMPLEX(SPC), DIMENSION(M), INTENT(OUT) :: roots 
 LOGICAL(LGT), INTENT(IN) :: polish 
 REAL(SP), PARAMETER :: EPS=1.0e-6_sp 
 !Given the array of M +1 complex coefficients a of the polynomial   
 !M+1 i=1 a i )x i  1 ,this routine successively calls laguer and  
 !finds all M complex roots The logical variable polish should be 
 !input as .true. if polishing (also by Laguerre's method) is 
 !desired, .false. if the roots will be subsequently polished by 
 !other means. Parameter: EPS is a small number. 
 INTEGER(I4B) :: j,its, I 
 COMPLEX(SPC) :: x 
 COMPLEX(SPC), DIMENSION(SIZE(a)) :: ad 
 !
  ad=a 
  !Copy of coefficients for successive deviation. 
  DO j=m,1,-1 
    !Loop over each root to be found. 
    x=CMPLX(0.0_sp,kind=spc)
    !Start at zero to favor convergence to smallest remaining root. 
    CALL laguer(ad(1:j+1),M,x,its) 
    !Find the root. 
    IF (ABS(AIMAG(x)) <= 2.0_sp*EPS**2*ABS(REAL(x))) & 
    x=CMPLX(REAL(x),kind=spc) 
    roots(j)=x 
    ad(j:1:-1)=poly_term_cc(ad(j+1:2:-1),x) 
    !Forward deflation. 
  END DO 
  IF (polish) THEN 
    DO j=1,m 
      !Polish the roots using the undeflated coefficients.
      CALL laguer(a(:),M,roots(j),its) 
    END DO 
  END IF 
   DO J=2,M
     X=ROOTS(J)
     DO I=J-1,1,-1
       IF(REAL(ROOTS(I)).LE.REAL(X))THEN
         ROOTS(I+1)=X
       ELSE
         ROOTS(I+1)=ROOTS(I)
       END IF
     END DO
     I=0
   END DO
!
 CONTAINS

   FUNCTION assert_eq2(n1,n2,string) 
!  ********************************* 
   !Report and die if integers not all equal (used for size checking). 
   CHARACTER(LEN=*), INTENT(IN) :: string 
   INTEGER, INTENT(IN) :: n1,n2 
   INTEGER :: assert_eq2 
    IF (n1 == n2) THEN 
      assert_eq2=n1 
    ELSE 
      WRITE (*,*)  'nrerror: an assert_eq failed with this tag:' , & 
        string 
      STOP  'program terminated by assert_eq2'  
    END IF
   END FUNCTION assert_eq2 

   RECURSIVE FUNCTION poly_term_cc(a,b) RESULT(u) 
   COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: a 
   COMPLEX(SPC), INTENT(IN) :: b 
   COMPLEX(SPC), DIMENSION(size(a)) :: u 
   INTEGER(I4B) :: n,j 
    n=size(a) 
    IF (n <= 0) RETURN 
    u(1)=a(1) 
    IF (n < 8) THEN 
      DO j=2,n 
        u(j)=a(j)+b*u(j-1) 
      END DO 
    ELSE 
      u(2:n:2)=poly_term_cc(a(2:n:2)+a(1:n-1:2)*b,b*b) 
      u(3:n:2)=a(3:n:2)+b*u(2:n-1:2) 
    END IF 
   END FUNCTION poly_term_cc
!
 END SUBROUTINE zroots
!____________________________________________________________________

 SUBROUTINE laguer(a,M,x,its) 
!****************************  
 USE OutPut, ONLY: Output_dir, file_ident
 IMPLICIT NONE 
 INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
 INTEGER, PARAMETER :: SP = KIND(1.0) 
 INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
 !
 INTEGER(I4B), INTENT(OUT) :: its 
 INTEGER(I4B), INTENT(IN) :: M
 COMPLEX(SPC), INTENT(INOUT) :: x 
 COMPLEX(SPC), DIMENSION(M+1), INTENT(IN) :: a 
 REAL(SP), PARAMETER :: EPS=epsilon(1.0_sp)*10 !*10;hack red. precision 
 INTEGER(I4B), PARAMETER :: MR=8,MT=10,MAXIT=MT*MR
 !Given an array of M +1 complex coe  cients a of the polynomial  
 !M+1 i=1 a i )x i  1 ,and given a complex value x this routine 
 !improves x by Laguerre's method until it converges, within 
 !the achievable roundoff limit, to a root of the given 
 !polynomial. The number of iterations taken is returned as its
 INTEGER(I4B) :: iter 
 REAL(SP) :: abx,abp,abm,err 
 COMPLEX(SPC) :: dx,x1,f,g,h,sq,gp,gm,g2 
 COMPLEX(SPC), DIMENSION(size(a)) :: b,d 
 REAL(SP), DIMENSION(MR) :: frac = &  !Fractions; to break limit cycle.
 (/ 0.5_sp,0.25_sp,0.75_sp,0.13_sp,0.38_sp,0.62_sp,0.88_sp,1.0_sp /)
 INTEGER, SAVE :: argh
  DO iter=1,MAXIT   !Loop over iterations up to allowed maximum.
    its=iter 
    abx=ABS(x) 
    b(m+1:1:-1)=poly_term_cc(a(m+1:1:-1),x) !compute  polynomial
    d(m:1:-1)=poly_term_cc(b(m+1:2:-1),x)   !and first two derivatives.
    f=poly_cc(x,d(2:m))                !f stores P /2 .
    err=EPS*poly_rr(abx,abs(b(1:m+1)))   !Estimate of roundoff.  
    IF (ABS(b(1)) <= err) RETURN   !We are on the root.
    g=d(1)/b(1)               !The generic case:Use Laguerre's formula.
    g2=g*g 
    h=g2-2.0_sp*f/b(1) 
    sq=SQRT((m-1)*(m*h-g2)) 
    gp=g+sq 
    gm=g-sq 
    abp=ABS(gp) 
    abm=ABS(gm) 
    IF (abp < abm) gp=gm 
    IF (MAX(abp,abm) > 0.0) THEN 
      dx=m/gp 
    ELSE 
      dx=EXP(CMPLX(LOG(1.0_sp+abx),iter,kind=spc)) 
    END IF 
    x1=x-dx
    IF (x == x1) RETURN   !Converged.
    IF (MOD(iter,MT) /= 0) THEN 
     x=x1
    ELSE                     !Every so often we take a fractional step, 
     x=x-dx*FRAC(iter/MT)    !to break any limit cycle (itself a rare 
    END IF                   !occurrence)
  END DO 
!Very unusual can occur only for complex roots.
!Try a different starting guess for the root.
  argh=argh+1
  OPEN(unit=21, file=TRIM(Output_dir)//'LAGUER_error_'//file_ident, & 
       status='replace')
    write(unit=21,fmt=*)argh, 'SUBROUTINE LAGUER: too many iterations in root find procedure'
  ENDFILE(unit=21) 
  CLOSE(unit=21) 
! 
!--------------------------------------------------------------------
CONTAINS
!  
  RECURSIVE FUNCTION poly_term_cc(a,b) RESULT(u)
! **********************************************   
  COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: a 
  COMPLEX(SPC), INTENT(IN) :: b 
  COMPLEX(SPC), DIMENSION(size(a)) :: u 
  INTEGER(I4B) :: n,j 
   n=size(a) 
   IF (n <= 0) RETURN 
   u(1)=a(1) 
   IF (n < 8) THEN 
     DO j=2,n 
       u(j)=a(j)+b*u(j-1) 
     END DO 
   ELSE 
     u(2:n:2)=poly_term_cc(a(2:n:2)+a(1:n-1:2)*b,b*b) 
     u(3:n:2)=a(3:n:2)+b*u(2:n-1:2) 
   END IF 
  END FUNCTION poly_term_cc
!--------------------------------------------------------------------
  FUNCTION poly_cc(x,coeffs)
! ************************** 
  COMPLEX(SPC), INTENT(IN) :: x 
  COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: coeffs 
  COMPLEX(SPC) :: poly_cc 
  COMPLEX(SPC) :: pow 
  COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: vec 
  INTEGER(I4B) :: i,n,nn 
   n=size(coeffs) 
   IF(n <= 0)THEN 
     poly_cc=0.0_sp 
   ELSE IF(n < 8)THEN 
     poly_cc=coeffs(n) 
     DO i=n-1,1,-1 
       poly_cc=x*poly_cc+coeffs(i) 
     END DO 
   ELSE 
     ALLOCATE(vec(n+1)) 
     pow=x 
     vec(1:n)=coeffs 
     DO 
       vec(n+1)=0.0_sp 
       nn=ishft(n+1,-1) 
       vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2) 
       IF(nn == 1)EXIT 
       pow=pow*pow 
       n=nn 
     END DO 
     poly_cc=vec(1) 
     DEALLOCATE(vec) 
   END IF 
  END FUNCTION poly_cc
!--------------------------------------------------------------------
  FUNCTION poly_rr(x,coeffs)
! **************************  
  !Polynomial evaluation. 
  REAL(SP), INTENT(IN) :: x 
  REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs 
  REAL(SP) :: poly_rr 
  REAL(SP) :: pow 
  REAL(SP), DIMENSION(:), ALLOCATABLE :: vec 
  INTEGER(I4B) :: i,n,nn 
   n=size(coeffs) 
   IF(n <= 0)THEN 
     poly_rr=0.0_sp 
   ELSE IF(n < 8)THEN 
     poly_rr=coeffs(n) 
     DO i=n-1,1,-1 
       poly_rr=x*poly_rr+coeffs(i) 
     END DO 
   ELSE 
     ALLOCATE(vec(n+1)) 
     pow=x 
     vec(1:n)=coeffs 
     DO 
       vec(n+1)=0.0_sp 
       nn=ishft(n+1,-1) 
       vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2) 
       IF(nn == 1)EXIT 
       pow=pow*pow 
       n=nn 
     END DO 
     poly_rr=vec(1) 
     DEALLOCATE(vec) 
   END IF 
  END FUNCTION poly_rr
!--------------------------------------------------------------------
 END SUBROUTINE laguer

  SUBROUTINE Flux
! ***************                                               
  USE Global, ONLY: L, c, pi
  USE Shock, ONLY: Radial_Height, b, RR, xR, R, x, Rrange, & 
                   xrange, rip_pck_num, SX, SY, SY_temp, RHSX
  USE Beam, ONLY: VR, vend, vcnt, vendold, vcntold, vmax,  & 
                  vcarray, vparray, f
  USE Spectra, ONLY: DS1F, DS2F, DS3F, DS1H, DS2H, DS3H, freq_bins,     & 
                     time_window, time_bins, exp_min, exp_range, l_val, & 
                     spec_shift, fber, loglem, blocked1, blocked2,      & 
                     blocked3, Xo1, Xo2, Xo3, Yo1, Yo2, Yo3
  USE Wind, ONLY: Vd, stheta1on, stheta, ctheta, Ne
  USE output, ONLY: total_F1, total_H1, total_F2, total_H2, total_F3, & 
                    total_H3
  IMPLICIT NONE 
! 
  REAL(KIND=L), ALLOCATABLE, DIMENSION(:) :: Yr, Zrsq
  REAL(KIND=L), ALLOCATABLE, DIMENSION(:) :: distsq1, distsq2, distsq3 
  INTEGER, ALLOCATABLE, DIMENSION(:) :: k1, k2, k3
  REAL(KIND=L) :: jH, jF 
  REAL(KIND=L) :: GRx, Rstep, xstep, Scale_2_3D, Areaofcalc
  REAL(KIND=L) :: vcoff
  REAL(KIND=L) :: jH_freq, jF_freq
  REAL(KIND=L) :: temp4aoc
  INTEGER :: Z, Q, i, j, rn, rpn, Alloc_error, Dealloc_error
! 
  rpn=rip_pck_num/2   ;   IF(rpn < 1)rpn=1
  Rstep=Rrange/(RR*RR) 
  xstep=xrange/(xR*xR)
  Scale_2_3D=1.2_L/b  !(0.3+0.3)/b=3Dising, 2.0=approx opposite theta 
                      !should incorporate explicit summation over slices 
                      !for detailed quantitative, however this is 
                      !generally OK to a factor ~few.
  i=fber*(LOG10(17.592_L*SQRT(Ne))-loglem) !2*fp(local plasma frequency)
!
  ALLOCATE(distsq1(rpn), distsq2(rpn), distsq3(rpn), STAT=alloc_error)
  IF(alloc_error /= 0)THEN 
    print*,'Failed to ALLOCATE arrays for distsq1, distsq2, distsq3'
    STOP
  END IF
!
  ALLOCATE(k1(rpn), k2(rpn), k3(rpn), STAT=alloc_error)
  IF(alloc_error /= 0)THEN 
    print*,'Failed to ALLOCATE arrays for k1, k2, k3'
    STOP
  END IF
!
  ALLOCATE(Yr(rpn), Zrsq(rpn), STAT=alloc_error)
  IF(alloc_error /= 0)THEN 
    print*,'Failed to ALLOCATE arrays for Yr, Zr'
    STOP
  END IF
!
  Yr(1)=SY_Temp
  DO rn=2,rpn
    Yr(rn)=SY_Temp*COS((rn-1)*(pi/2)/(rpn-1)) 
  END DO
  Zrsq=SY*SY - Yr*Yr
!
  distsq1=(Xo1-RHSX)**2 + (Yo1-Yr)**2 + Zrsq
  distsq2=(Xo2-RHSX)**2 + (Yo2-Yr)**2 + Zrsq
  distsq3=(Xo3-RHSX)**2 + (Yo3-Yr)**2 + Zrsq
! convert distance from observer into a time bin number:
  k1=1+(SQRT(distsq1)/c)/((1.0_L*time_window)/time_bins) 
  k2=1+(SQRT(distsq2)/c)/((1.0_L*time_window)/time_bins) 
  k3=1+(SQRT(distsq3)/c)/((1.0_L*time_window)/time_bins) 
! time bin based on straight line
! Above assumes ripple is far enough away that geometry of source 
! region is unimportant & thus total volume emissivity divided by 
! (one) distance squared will give good approx to flux...
!
  vend=VR
  vcnt=1
!
! Do for array of points over a 2D shock/foreshock region:  
!
  DO Z=1,RR
    R=Z*Z*Rstep
    temp4aoc=xstep*(2.0_L*Z-1.0_L)*Rstep
!   counting backward in x position so that previous Vc
!   may be used as starting point for following calculation.
    DO Q=xR,1,-1
     x=Q*Q*xstep 
     Areaofcalc=(2.0_L*Q-1.0_L)*temp4aoc
!
!    Retain previous range of beam for limiting calculation range:
!
     vcntold=MAX(1,vcnt-3) ; vendold=vend
!
!    Initialise cutoff velocity counter to 1 for first calculation
!    in each column.   
!
     IF(Q == xR)vcnt=1 
!  
!    Cut out of calculations regions where velocity is to high 
!    or to low such that it will cause problems in the qlf subroutine.
!    This is a physically acceptable removal as it corresponds to
!    velocities at which physical restraints would remove them anyway.
!
     vcoff=Vd*R/x
     IF(vcoff < vmax)THEN
!     
!    Only calculate for positions upstream of shock:
! 
       GRx=x*stheta1on+b*(x*ctheta-R*stheta)**2
       IF(GRx > 0.0_L)THEN
         CALL cutoff_velocity_distribution
!
!        Combine the background and cutoff distributions
!        and calculate the quasi-linear flattening effect:
!
         f(vcntold:vendold)=vparray(vcntold:vendold)
         f(vcnt:vend)=f(vcnt:vend)+vcarray(vcnt:vend)
         CALL qlf
!
!        Calculate the volume emissivity from this location:
!
         CALL emissivities (jH,jH_freq,jF,jF_freq)
!
       ELSE
!      PRINT*, 'Inside Shock'
         jH=0.0_L ; jF=0.0_L
       ENDIF  
     ELSE
!    PRINT*, 'cutoff velocity too high or low for useful beam to form'
       jH=0.0_L ; jF=0.0_L
     ENDIF
!  
! IF blocked use frequency restriction requirement: 
!
     IF(jH > 0.0)THEN
      jH=jH*Areaofcalc*Scale_2_3D
      j=fber*(LOG10(jH_freq)-loglem)  !index for Frequency of emission
      IF(j > freq_bins)j=freq_bins
      IF(j < 1)j=1
!
      IF(blocked1 == 0)THEN
       DO rn=1,rpn
        DS1H(j,k1(rn):k1(rn)+spec_shift)=DS1H(j,k1(rn):k1(rn)+spec_shift) & 
                                         +2*jH/distsq1(rn)
       END DO
       total_H1=total_H1+2*jH/distsq1(1)
      ELSE 
       IF(j > i)THEN 
        DO rn=1,rpn
         DS1H(j,k1(rn):k1(rn)+spec_shift)=DS1H(j,k1(rn):k1(rn)+spec_shift) & 
                                          +2*jH/distsq1(rn)
        END DO
        total_H1=total_H1+2*jH/distsq1(1)
       END IF
      END IF
!
      IF(blocked2 == 0)THEN
       DO rn=1,rpn
        DS2H(j,k2(rn):k2(rn)+spec_shift)=DS2H(j,k2(rn):k2(rn)+spec_shift) & 
                                         +2*jH/distsq2(rn)
       END DO
       total_H2=total_H2+2*jH/distsq2(1)
      ELSE 
       IF(j > i)THEN 
        DO rn=1,rpn
         DS2H(j,k2(rn):k2(rn)+spec_shift)=DS2H(j,k2(rn):k2(rn)+spec_shift) & 
                                          +2*jH/distsq2(rn)
        END DO
        total_H2=total_H2+2*jH/distsq2(1)
       END IF
      END IF
!
      IF(blocked3 == 0)THEN
       DO rn=1,rpn
         DS3H(j,k3(rn):k3(rn)+spec_shift)=DS3H(j,k3(rn):k3(rn)+spec_shift) & 
                                          +2*jH/distsq3(rn)
       END DO
       total_H3=total_H3+2*jH/distsq3(1)
      ELSE 
       IF(j > i)THEN 
        DO rn=1,rpn
         DS3H(j,k3(rn):k3(rn)+spec_shift)=DS3H(j,k3(rn):k3(rn)+spec_shift) & 
                                          +2*jH/distsq3(rn)
        END DO
        total_H3=total_H3+2*jH/distsq3(1)
       END IF
      END IF
     END IF
!
     IF(jF > 0.0)THEN
      jF=jF*Areaofcalc*Scale_2_3D
      j=fber*(LOG10(jF_freq)-loglem) !index for Frequency of emission
      IF(j > freq_bins)j=freq_bins
      IF(j < 1)j=1
!
      IF(blocked1 == 0)THEN
       DO rn=1,rpn
        DS1F(j,k1(rn):k1(rn)+spec_shift)=DS1F(j,k1(rn):k1(rn)+spec_shift) & 
                                         +2*jF/distsq1(rn)
       END DO
       total_F1=total_F1+2*jF/distsq1(1)
      ELSE 
       IF(j > i)THEN 
        DO rn=1,rpn
         DS1F(j,k1(rn):k1(rn)+spec_shift)=DS1F(j,k1(rn):k1(rn)+spec_shift) & 
                                          +2*jF/distsq1(rn)
        END DO
        total_F1=total_F1+2*jF/distsq1(1)
       END IF
      END IF
!
      IF(blocked2 == 0)THEN
       DO rn=1,rpn
        DS2F(j,k2(rn):k2(rn)+spec_shift)=DS2F(j,k2(rn):k2(rn)+spec_shift) & 
                                         +2*jF/distsq2(rn)
       END DO
       total_F2=total_F2+2*jF/distsq2(1)
      ELSE 
       IF(j > i)THEN 
        DO rn=1,rpn
         DS2F(j,k2(rn):k2(rn)+spec_shift)=DS2F(j,k2(rn):k2(rn)+spec_shift) & 
                                          +2*jF/distsq2(rn)
        END DO
        total_F2=total_F2+2*jF/distsq2(1)
       END IF
      END IF
!
      IF(blocked3 == 0)THEN
       DO rn=1,rpn
        DS3F(j,k3(rn):k3(rn)+spec_shift)=DS3F(j,k3(rn):k3(rn)+spec_shift) & 
                                         +2*jF/distsq3(rn)
       END DO
       total_F3=total_F3+2*jF/distsq3(1)
      ELSE 
       IF(j > i)THEN 
        DO rn=1,rpn
         DS3F(j,k3(rn):k3(rn)+spec_shift)=DS3F(j,k3(rn):k3(rn)+spec_shift) & 
                                          +2*jF/distsq3(rn)
        END DO
        total_F3=total_F3+2*jF/distsq3(1)
       END IF
      END IF         
     END IF
!
    END DO
  END DO
!  
  DEALLOCATE(distsq1, distsq2, distsq3, STAT=dealloc_error)
  IF(dealloc_error /= 0)PRINT*,'Deallocation error for distsq arrays'
!  
  DEALLOCATE(k1, k2, k3, STAT=dealloc_error)
  IF(dealloc_error /= 0)PRINT*,'Deallocation error for k arrays'
!  
  DEALLOCATE(Yr, Zrsq, STAT=dealloc_error)
  IF(dealloc_error /= 0)PRINT*,'Deallocation error for Yr, Zrsq arrays'
!
  END SUBROUTINE Flux 
!____________________________________________________________________

  SUBROUTINE cutoff_velocity_distribution
! ***************************************                      
  USE Global, ONLY: c, csq
  USE Beam
  USE Wind, ONLY: stheta, stheta2, ctheta, Vesq, kappa
  USE Shock, ONLY: R, x, SRonSrange, Vc, alph, B2B1, B2B1_1, &
                   UnSalphonB2B1, ephi2onme, S_limit
  IMPLICIT NONE
!
! Rs= R coordinate on shock from which electron with particular velocity
! would have come; xs=x as for R; Vpi=initial parallel velocity;
! vperlcsq= square of perp velocity below which electrons cannot escape 
! from shock based on losscone; vpersdsq= as above but based on shock
! drift acceleration constraints; vpersq= dominant one of previous two;  
!
  REAL(KIND=L) :: Rs, xs, Vpi, vperlcsq, vpersdsq, vpersq, A, &
                  Vpisq, shockx, temp4vcarr, vcmax, Atemp4Rvponvdx
!   
  INTEGER :: cntr, S, i
!
  vcarray(vcntold:vendold)=0.0_L
  vcmax=0.0_L
  cntr=0
!
  DO i=vcnt,VR
    A=ctheta-temp4vponvd(i)*stheta
    Atemp4Rvponvdx=A*(R-temp4vponvd(i)*x)
    shockx=bsqsthetasq-(b4on*Atemp4Rvponvdx)
!    
!   Elimination of distribution due to complex results in travel
!   and shockx are equivalent to a cut off velocity due to time
!   of flight:
!   
    IF(shockx < 0.0_L)CYCLE 
! 
!   Determine position on shock from which a particle with a
!   particular parallel velocity has come:
!     
    xs=(stheta2*Atemp4Rvponvdx-bstheta1on+SQRT(shockx))/(2.0*A*A)
!
!   Require that the position is upstream of current location
!   and therefore physical:
!
    IF(xs > x)CYCLE 
!
    Rs=R-temp4vponvd(i)*(x-xs)
    S=SRonSrange*Rs
    IF(S >= S_limit)THEN 
      vend=i
      RETURN !related to thetaBn/alph
    END IF 
!         
!   Initiate counter to limit do loops (ie prevent calculations
!   below vc): 
!
    IF(cntr == 0)THEN 
      vcnt=Vc(S)/dv
      cntr=1
    END IF
! 
    IF(alph(S) <= 1.0e-13)CYCLE  !Only really need for high resolution 
                                 !where the shock is sampled near the 
                                 !tangent point.
!
!   Calculate initial parallel speed from which reflected electrons 
!   came:
!
    Vpi=Uctp+2.0_L*Vc(S)-vel(i)
!   
!   Require that incident electrons encounter the shock and 
!   travel at < the speed of light: 
!
    IF((Vpi > Vc(S)) .OR. (Vpi < -c))CYCLE 
!   
!   Incorporate requirement for reflection via shock drift
!   acceleration to occur:
!
    IF(B2B1(S) > 1.0_L)THEN 
      Vpisq=Vpi*Vpi 
      vperlcsq=((Vpi-Vc(S))**2+ephi2onme(S))/B2B1_1(S)
      vpersdsq=UnSalphonB2B1(S)-Vpisq        
      vpersq=MAX(vperlcsq,vpersdsq)
!         
!     Calculate reflected distribution:
!     
      IF(vpersq < csq)THEN  
        temp4vcarr=Vpisq+Vesq
        vcarray(i)=ConstDist*(((temp4vcarr+vpersq)*1.0e-6)**(-Kappa)  &
                   -((temp4vcarr+csq)*1.0e-6)**(-Kappa))
      END IF
    END IF
!   
!   Stop calculating when number of reflected electrons drops 
!   below 1% of the peak height of the reflected distribution
!   (valid when only one beam is present, ie smooth reflection 
!   and acceleration at shock):
!     
    IF(vcarray(i) > vcmax)THEN
      vcmax=vcarray(i)
    ELSE IF(vcarray(i) < 0.01_L*vcmax)THEN 
      vend=i
      RETURN
    END IF
!
  END DO
!
  END SUBROUTINE cutoff_velocity_distribution
!____________________________________________________________________

  SUBROUTINE qlf
! ***************
!
! Estimates beam parameters.
!                                               
  USE Global, ONLY: L
  USE Beam, ONLY: deltavb, vb, Nb, vcnt, VR, f, dv, halfdv, vcarray
  USE Wind, ONLY: Ne
  IMPLICIT NONE
  REAL(KIND=L) :: dfdv, area, d1, d2
  INTEGER :: w, p1, p2, p3
!
! Initialize beam parameters:
!      
  deltavb=0.0_L  ! flattened beam width
  vb=0.0_L       ! flattened beam speed
  Nb=0.0_L
  area=0.0_L
!
! Calculate flattened beam distribution:
!  
  p2=MAX(vcnt-3,0)  
!
  DO 
    p2=p2+1
    IF(p2 >= VR-3)EXIT
    dfdv=f(p2+1)-f(p2)
    IF(dfdv > 0.0)THEN
      w=1 ; p1=p2 ; p3=p2+1 ; p2=p2+2 
      area=f(p3) 
      d1=area-f(p1) ; d2=f(p2)-area 
      DO 
        IF((d1 > d2).AND.(d1 > 0.0))THEN 
          w=w+1 ; area=area+f(p1) ; p1=p1-1
      IF(p1 < 1) EXIT
          d1=area-f(p1)*w ; d2=f(p2)*w-area
        ELSE IF(d2 > 0.0)THEN
          w=w+1 ; area=area+f(p2) ; p2=p2+1 
      IF(p2 > VR) EXIT
          d1=area-f(p1)*w ; d2=f(p2)*w-area
        ELSE 
          EXIT 
        END IF 
      END DO
      vb=halfdv*(p2+p1)
      deltavb=halfdv*(p2-p1-2)
      Nb=Ne*dv*SUM(vcarray(p1+1:p2-1))
      EXIT                             
    END IF                            
  END DO
!
  END SUBROUTINE qlf
!____________________________________________________________________

  SUBROUTINE emissivities (jH1,H_freq,jF1,F_freq)
! *************************************************** 
  USE Global, ONLY: L, sqtwo, c, sqpis, me
  USE Emis, ONLY: temp4gammaLS, temp4zetaF, Vcrit, delohmF, delohmH, & 
                  erftemp1, erftemp2, temp4phiF, temp4phiH, erf,     & 
                  gmemisq, Gaussian_Random_e
  USE Shock, ONLY: R,x
  USE Beam, ONLY: deltavb, vb, Nb
  USE Wind, ONLY: Vesq, Ve, Ne
  IMPLICIT NONE
!
! See Robinson & Cairns 1998, Solar Physics 181, 395.
!
  REAL(KIND=L), INTENT(OUT):: jH1,jF1
  REAL(KIND=L), INTENT(OUT):: H_freq,F_freq
  REAL(KIND=L), DIMENSION(1) :: BOB=0.0_L
  REAL(KIND=L), PARAMETER:: alphas=1.0_L, beta=0.3333_L
  REAL(KIND=L):: phiF, phiH, zetaF, zetaH, gammaLS, &
                 lr, temp4j, temp4zetaH, vbsq, fp
!
  IF((deltavb == 0.0_L).OR.(vb == 0.0_L))THEN
    jF1=0.0_L  ;  jH1=0.0_L
  ELSE
    vbsq=vb*vb
    lr=SQRT(R*R+x*x)
    gammaLS=temp4gammaLS/(alphas*vbsq)
    zetaF=EXP(-temp4zetaF*(vb/(beta*deltavb))**2 &
          *((Vcrit-vb)**2/Vesq))
    temp4zetaH=Ve*beta*deltavb*sqtwo/vbsq
    zetaH=(c*sqpis*beta*deltavb/(2.0_L*vbsq)) &
          *(erf(erftemp1/temp4zetaH)          &
          +erf(erftemp2/temp4zetaH))
    phiF=temp4phiF*vb*zetaF*gammaLS
    phiH=temp4phiH*zetaH*vbsq*vb
    temp4j=me*Nb*vbsq/(3.0_L*lr)
    jF1=(phiF/delohmF)*temp4j
    jH1=(phiH/delohmH)*temp4j  
    IF(vb < 0.7_L*Ve)jF1=0.0_L 
!
! If emission is produced calculate its frequency using a 
! Gaussian random fluctuation in Ne: 
! 
    IF((jF1 > 0.0) .OR. (jH1 > 0.0))THEN
      CALL Gaussian_Random_e(BOB) 
      IF(bob(1) < -1.0)bob(1)=-1.0_L
      IF(bob(1) > 1.0)bob(1)=1.0_L 
      fp=8.976_L*SQRT((1.0_L+(0.2_L*bob(1)))*Ne) 
    END IF 
    IF(jF1 > 0.0)F_freq=fp*(1.0_L+(1.5_L*Vesq/vbsq)-(Ve/vb)*gmemisq)
    IF(jH1 > 0.0)H_freq=2.0_L*fp 
!    
  END IF
!    
  END SUBROUTINE emissivities
!____________________________________________________________________

  SUBROUTINE Finalise_and_Write
! *****************************       
  USE Global, ONLY: L
  USE Spectra, ONLY: DS1F, DS2F, DS3F, DS1H, DS2H, DS3H, & 
                     exp_min, freq_bins, spec_shift, erfb
  USE Shock, ONLY: UT_time 
  IMPLICIT NONE
!  
  REAL(KIND=L) :: divby, temp4divby
  INTEGER :: i, j
!
! divides by frequency width of bin to get flux per unit frequency
! 
  temp4divby=1.0_L-10.0_L**(-erfb)
!
  DO i=1,freq_bins 
    divby=10.0_L**(exp_min+i*erfb)*temp4divby
!
    DO j=1,spec_shift 
      IF(DS1F(i,j) > 1.e-38_L)THEN 
        DS1F(i,j)=LOG10(DS1F(i,j)/divby) 
      ELSE 
        DS1F(i,j)=0.0_L 
      END IF   
      IF(DS1H(i,j) > 1.e-38_L)THEN 
        DS1H(i,j)=LOG10(DS1H(i,j)/divby) 
      ELSE 
        DS1H(i,j)=0.0_L 
      END IF  
    END DO 
    DO j=1,spec_shift 
      IF(DS2F(i,j) > 1.e-38_L)THEN 
        DS2F(i,j)=LOG10(DS2F(i,j)/divby) 
      ELSE 
        DS2F(i,j)=0.0_L 
      END IF  
      IF(DS2H(i,j) > 1.e-38_L)THEN 
        DS2H(i,j)=LOG10(DS2H(i,j)/divby) 
      ELSE 
        DS2H(i,j)=0.0_L 
      END IF  
    END DO 
    DO j=1,spec_shift 
      IF(DS3F(i,j) > 1.e-38_L)THEN 
        DS3F(i,j)=LOG10(DS3F(i,j)/divby) 
      ELSE 
        DS3F(i,j)=0.0_L 
      END IF  
      IF(DS3H(i,j) > 1.e-38_L)THEN 
        DS3H(i,j)=LOG10(DS3H(i,j)/divby) 
      ELSE 
        DS3H(i,j)=0.0_L 
      END IF  
    END DO 
!     
  END DO
!
  CALL Write_to_file 
!
! Thow away the stuff just written to a file:
!
  DS1F=EOSHIFT(DS1F,SHIFT=spec_shift,BOUNDARY=0.0_L,DIM=2)
  DS2F=EOSHIFT(DS2F,SHIFT=spec_shift,BOUNDARY=0.0_L,DIM=2)
  DS3F=EOSHIFT(DS3F,SHIFT=spec_shift,BOUNDARY=0.0_L,DIM=2) 
  DS1H=EOSHIFT(DS1H,SHIFT=spec_shift,BOUNDARY=0.0_L,DIM=2)
  DS2H=EOSHIFT(DS2H,SHIFT=spec_shift,BOUNDARY=0.0_L,DIM=2)
  DS3H=EOSHIFT(DS3H,SHIFT=spec_shift,BOUNDARY=0.0_L,DIM=2)
!
  END SUBROUTINE Finalise_and_Write
!____________________________________________________________________

  SUBROUTINE Write_to_file
! ************************
  USE OutPut
  USE Wind, ONLY: Average_kappa, Vsw, nos, SWS, Ne_1AU, Ne_1Rs, Te_1AU, & 
                  Ti_1AU, Bo, Vsw_1AU
  USE Shock, ONLY: Radial_Height, UT_time, UCME, CME_accel, & 
                   expansion, scale, scale_exp
  USE Spectra, ONLY: DS1F, DS2F, DS3F, DS1H, DS2H, DS3H,    & 
                     exp_min, exp_range, l_val, spec_shift, & 
                     Xo1, Xo2, Xo3, Yo1, Yo2, Yo3
  IMPLICIT NONE 
! 
  CHARACTER(LEN=7) :: output_status
  CHARACTER(LEN=77) :: full_path_0a, full_path_0b, full_path_0c, & 
                       full_path_0d, full_path_0e, full_path_0f, & 
                       full_path_1, full_path_1b, full_path_2,   & 
                       full_path_3, full_path_4, full_path_5,    &
                       full_path_6, full_path_8, full_path_9,    &
                       full_path_10, full_path_11
  INTEGER :: j
!    
  330 FORMAT('********************************',/,                      &
             '*** INITIAL PARAMETER VALUES ***',/,                      &
             '********************************',//,                     &
             'Radial_Height=',EN14.5,' m ',//,                          &
             'CME Speed=',EN13.4,' m/s',//,                             &
             'CME deceleration=',F4.1,' m/s/s ',//,                     &
             'CME expansion=',F4.1,                                     & 
                            ' (0=no expansion;1=expanding at UCME)',//, &  
             'CME scale=',F4.1,' (1.0=1AU; Bigger=smaller)',//,         &  
             'CME scale rate of change=',F4.2,                          & 
               ' (for scale=1~>0.30; scale=3~>0.37)',//,                & 
             'kappa=',F4.1,//,                                          &
             '1AU normalizations:',/,                                   &
             'Ne=',EN11.3,' m^-3',/,                                    &
             'Ne=',EN11.3,' m^-3  (coronal value)',/,                   &
             'Te=',EN11.3,' K',/,                                       &
             'Ti=',EN11.3,' K',/,                                       &
             'Bo=',EN16.8,' T  (coronal value)',/,                      &
             'Vsw=',EN11.3,' m/s',//,                                   &
             'Exp_min=',F4.1,/,                                         &
             'Exp_range=',F4.1,/,                                       &
             'l_val=',F4.1, ///,                                        &
             'Observer 1',/,                                            &
             'X=',ES11.2,' m',/,                                        &
             'Y=',ES12.2,' m',//,                                       &
             'Observer 2',/,                                            &
             'X=',ES11.2,' m',/,                                        &
             'Y=',ES12.2,' m',//,                                       &
             'Observer 3',/,                                            &
             'X=',ES11.2,' m',/,                                        &
             'Y=',ES12.2,' m',/// )   
! 
  370 FORMAT('*********************************',/,   & 
             '*** EXPLICIT STRUCTURE VALUES ***',/,   & 
             '*********************************',//,  & 
             'There are ',I2,' explicitly defined',/, & 
             'structures in this calculation.',/,     & 
             'They are:',// )
  371 FORMAT('Structure ',I2,' is of type ', I2,': a CIR',/,        & 
             'Rotation=',F4.2,' ; ','Width=',F4.2,// ) 
  372 FORMAT('Structure ',I2,' is of type ', I2, & 
             ': a heliocentric annulus',/,       &
             'with inner and outer radius',/,    &
             'R_inner=',ES12.3,' m',/,           &
             'R_outer=',ES12.3,' m',/,           &
             'respectively.',// )  
  374 FORMAT('Structure ',I2,' is of type ', I2,  & 
             ': a square/rectangular',/,          &  
             'region bound by',/,                 &
             'X_min=',ES12.3,' m',' ; ','X_max=', & 
             ES12.3,' m',/,                       &
             'Y_min=',ES12.3,' m',' ; ','Y_max=', & 
             ES12.3,' m',// )  
  373 FORMAT('Structure ',I2,' is of type ', I2,            & 
             ': a circular region of radius',/,             &
             'R=',ES12.3,' m',/,                            &
             'and centered on heliocentric coordinates',/,  &
             'X=',ES12.3,' m',' ; ','Y=',ES12.3,' m',// )  
  376 FORMAT('Structure ',I2,' is of type ', I2,            &
             ': a coronal loop, defined',/,                 &
             'by two overlapping circles',/,                &
             'R1=',ES12.3,' m',/,                           &
             'X1=',ES12.3,' m',' ; ','Y1=',ES12.3,' m',/,   &
             'R2=',ES12.3,' m',/,                           &
             'X2=',ES12.3,' m',' ; ','Y2=',ES12.3,' m',// )
! 
!
  350 FORMAT(E13.7) 
!  
  300 FORMAT(/,'UT =',I7,' sec',/,'Vsw=',EN14.5,' m/s ')
!
  full_path_0a=TRIM(Output_dir)//'DS1F_'//file_ident
  full_path_0b=TRIM(Output_dir)//'DS2F_'//file_ident
  full_path_0c=TRIM(Output_dir)//'DS3F_'//file_ident
  full_path_0d=TRIM(Output_dir)//'DS1H_'//file_ident
  full_path_0e=TRIM(Output_dir)//'DS2H_'//file_ident
  full_path_0f=TRIM(Output_dir)//'DS3H_'//file_ident
  full_path_1=TRIM(Output_dir)//'parameters_'//file_ident
  full_path_1b=TRIM(Output_dir)//'scaling_'//file_ident
  full_path_2=TRIM(Output_dir)//'Radial_Height_'//file_ident
  full_path_3=TRIM(Output_dir)//'UCME_'//file_ident
  full_path_4=TRIM(Output_dir)//'theta_'//file_ident
  full_path_5=TRIM(Output_dir)//'U_'//file_ident
  full_path_6=TRIM(Output_dir)//'b_'//file_ident
  full_path_8=TRIM(Output_dir)//'B1_'//file_ident
  full_path_9=TRIM(Output_dir)//'Ne_'//file_ident
  full_path_10=TRIM(Output_dir)//'Te_'//file_ident
  full_path_11=TRIM(Output_dir)//'Ti_'//file_ident
! 
  IF(first_pass == 1)THEN 
    IF(file_ident == 'temp.dat')THEN
    
      OPEN(unit=21, file=TRIM(Output_dir)//'LOG_'//file_ident, & 
           status='old', position='append')
        write(unit=21,fmt=*) 'Writing over existing *temp.dat files...' 
      ENDFILE(unit=21) 
      CLOSE(unit=21)
      output_status='replace'
    ELSE
      output_status='new'
    END IF
!
    OPEN(unit=14,file=full_path_1,status=output_status)
      write(unit=14,fmt=330)Radial_Height, UCME, CME_accel, expansion, & 
                            scale, scale_exp, Average_kappa, Ne_1AU,   & 
                            Ne_1Rs, Te_1AU, Ti_1AU, Bo, Vsw_1AU,       & 
                            exp_min, exp_range, l_val,                 & 
                            Xo1, Yo1, Xo2, Yo2, Xo3, Yo3
      write(unit=14,fmt=370)nos 
      DO j=1,nos
        SELECT CASE (SWS(j)%shape)
          CASE (1)
            write(unit=14,fmt=371)j,SWS(j)%shape,SWS(j)%param(1), & 
                                  SWS(j)%param(2)
          CASE (2)
            write(unit=14,fmt=372)j,SWS(j)%shape,SWS(j)%param(1), & 
                                  SWS(j)%param(2)
          CASE (4)
            write(unit=14,fmt=374)j,SWS(j)%shape,SWS(j)%param(1),  & 
                                  SWS(j)%param(2),SWS(j)%param(3), & 
                                  SWS(j)%param(4)
          CASE (3)
            write(unit=14,fmt=373)j,SWS(j)%shape,SWS(j)%param(1), & 
                                  SWS(j)%param(2),SWS(j)%param(3)
          CASE (6)
            write(unit=14,fmt=376)j,SWS(j)%shape,SWS(j)%param(1),  & 
                                  SWS(j)%param(2),SWS(j)%param(3), & 
                                  SWS(j)%param(3),SWS(j)%param(4), & 
                                  SWS(j)%param(5),SWS(j)%param(6)
        END SELECT
      END DO 
    ENDFILE(unit=14) 
    CLOSE(unit=14)
!   
    OPEN(unit=13,file=full_path_0a,status=output_status) 
    OPEN(unit=15,file=full_path_0b,status=output_status) 
    OPEN(unit=16,file=full_path_0c,status=output_status)  
    ENDFILE(unit=13) 
    CLOSE(unit=13)
    ENDFILE(unit=15) 
    CLOSE(unit=15)
    ENDFILE(unit=16) 
    CLOSE(unit=16)
!   
    OPEN(unit=43,file=full_path_0d,status=output_status) 
    OPEN(unit=45,file=full_path_0e,status=output_status) 
    OPEN(unit=46,file=full_path_0f,status=output_status)  
    ENDFILE(unit=43) 
    CLOSE(unit=43)
    ENDFILE(unit=45) 
    CLOSE(unit=45)
    ENDFILE(unit=46) 
    CLOSE(unit=46)
!
    OPEN(unit=13,file=full_path_1b,status=output_status)
    ENDFILE(unit=13) 
    CLOSE(unit=13)
!    
    OPEN(unit=13,file=full_path_2,status=output_status)
    ENDFILE(unit=13) 
    CLOSE(unit=13)
!
    OPEN(unit=13,file=full_path_3,status=output_status)
    ENDFILE(unit=13) 
    CLOSE(unit=13)
!
    OPEN(unit=13,file=full_path_4,status=output_status)
    ENDFILE(unit=13) 
    CLOSE(unit=13)
! 
    OPEN(unit=13,file=full_path_5,status=output_status)
    ENDFILE(unit=13) 
    CLOSE(unit=13)
!
    OPEN(unit=13,file=full_path_6,status=output_status)
    ENDFILE(unit=13) 
    CLOSE(unit=13)
!
    OPEN(unit=13,file=full_path_8,status=output_status)
    ENDFILE(unit=13) 
    CLOSE(unit=13)
!
    OPEN(unit=13,file=full_path_9,status=output_status)
    ENDFILE(unit=13) 
    CLOSE(unit=13)
! 
    OPEN(unit=13,file=full_path_10,status=output_status)
    ENDFILE(unit=13) 
    CLOSE(unit=13) 
! 
    OPEN(unit=13,file=full_path_11,status=output_status)
    ENDFILE(unit=13) 
    CLOSE(unit=13)
!
    first_pass=0
!    
  END IF
! 
  OPEN(unit=13,file=full_path_0a,status='old',position='append')
  OPEN(unit=15,file=full_path_0b,status='old',position='append')
  OPEN(unit=16,file=full_path_0c,status='old',position='append')
  OPEN(unit=43,file=full_path_0d,status='old',position='append')
  OPEN(unit=45,file=full_path_0e,status='old',position='append')
  OPEN(unit=46,file=full_path_0f,status='old',position='append')
  OPEN(unit=14,file=full_path_1,status='old',position='append')
    write(unit=13,fmt=350)DS1F(:,1:spec_shift)
    write(unit=15,fmt=350)DS2F(:,1:spec_shift)
    write(unit=16,fmt=350)DS3F(:,1:spec_shift)
    write(unit=43,fmt=350)DS1H(:,1:spec_shift)
    write(unit=45,fmt=350)DS2H(:,1:spec_shift)
    write(unit=46,fmt=350)DS3H(:,1:spec_shift)
    write(unit=14,fmt=300)UT_time,Vsw
  ENDFILE(unit=13) 
  CLOSE(unit=13)
  ENDFILE(unit=15) 
  CLOSE(unit=15)
  ENDFILE(unit=16) 
  CLOSE(unit=16)
  ENDFILE(unit=43) 
  CLOSE(unit=43)
  ENDFILE(unit=45) 
  CLOSE(unit=45)
  ENDFILE(unit=46) 
  CLOSE(unit=46)
  ENDFILE(unit=14) 
  CLOSE(unit=14)
!
  OPEN(unit=13,file=full_path_1b,status='old',position='append')
    write(unit=13,fmt=*)spec_shift
  ENDFILE(unit=13) 
  CLOSE(unit=13)
!
  OPEN(unit=13,file=full_path_2,status='old',position='append')
  OPEN(unit=14,file=full_path_3,status='old',position='append')
    write(unit=13,fmt=*)Radial_Height
    write(unit=14,fmt=*)UCME
  ENDFILE(unit=13) 
  CLOSE(unit=13)
  ENDFILE(unit=14) 
  CLOSE(unit=14)
!
  OPEN(unit=13,file=full_path_4,status='old',position='append')
    write(unit=13,fmt=*)theta_arr
  ENDFILE(unit=13) 
  CLOSE(unit=13)
!
  OPEN(unit=13,file=full_path_5,status='old',position='append')
    write(unit=13,fmt=*)U_arr
  ENDFILE(unit=13) 
  CLOSE(unit=13)
!
  OPEN(unit=13,file=full_path_6,status='old',position='append')
    write(unit=13,fmt=*)b_arr
  ENDFILE(unit=13) 
  CLOSE(unit=13)
!
  OPEN(unit=13,file=full_path_8,status='old',position='append')
    write(unit=13,fmt=*)B1_arr
  ENDFILE(unit=13) 
  CLOSE(unit=14)
!
  OPEN(unit=13,file=full_path_9,status='old',position='append')
    write(unit=13,fmt=*)Ne_arr
  ENDFILE(unit=13) 
  CLOSE(unit=13)
!
  OPEN(unit=13,file=full_path_10,status='old',position='append')
    write(unit=13,fmt=*)Te_arr
  ENDFILE(unit=13) 
  CLOSE(unit=13)
!
  OPEN(unit=13,file=full_path_11,status='old',position='append')
    write(unit=13,fmt=*)Ti_arr
  ENDFILE(unit=13) 
  CLOSE(unit=13)
!  
  END SUBROUTINE Write_to_file          
!____________________________________________________________________
!
!EOF Dynamic_Spectra.f90  (Thesis Version)
