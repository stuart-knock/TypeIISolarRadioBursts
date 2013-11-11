# Type II Solar Radio Bursts

This code calculates the predicted dynamic spectra for a type II  radio burst
from the solar corona and interplanetary medium. The program is written in 
Fortran 90, formatted output is written to file with plotting of that output
performed by two IDL routines. 

It's the code from my PhD.

-------------------------------------------------------------------------------

## Quick Start Guide:

1) Change into the FORTRAN/ directory and compile Dynamic_Spectra.f90
    
    #Compiling with gcc to run code:
    
      gfortran -O3 -o Dynamic_Spectra.e Dynamic_Spectra.f90
    
    #Compiling with gcc to debug code:
    
      gfortran -Wall -Wextra -Wconversion -pedantic -g -o Dynamic_Spectra.e Dynamic_Spectra.f90
    
    #Compiling with the intel fortran compiler to run code:
    
      ifort -O3 -xHost -o Dynamic_Spectra.e Dynamic_Spectra.f90
    
    #Compiling with intel fortran compiler to debug code:

      ifort -debug full -debug emit_column -debug inline-debug-info -debug extended -check -traceback -fpe0 -o Dynamic_Spectra.e Dynamic_Spectra.f90


2) Run the code using the example input file by typing 
     
    Dynamic_Spectra.e

at a command line. The program uses input from Dynamic_Spectra.input if it can 
find it otherwise it uses default values hard coded into the subroutine 
Initialise_Global. The example Dynamic_Spectra.input file takes about 3 
minutes on a circa 2012 workstation  and produces about 42MB of .dat files. 

*Back when I first wrote this, the same default simulation took 30 minutes on 
a ~2GHz P4/Athlon.*

Once the code has finished running, run the Interactive Data Language (IDL)
routines in the IDL/ directory to plot the output.

-------------------------------------------------------------------------------

## Contents of this Directory:

**doc/**

Contains a LateX document which briefly describes the code.

 --------------------------------------------------------------------
  
**OUTPUT/**

The directory to which output from Dynamic_Spectra.f90 is written.

 --------------------------------------------------------------------
  
**FORTRAN/**

This directory contains the Fortran source code and an example of the input 
file it requires. 

 --------------------------------------------------------------------
  
**IDL/**

This directory contains IDL routines for plotting the output from 
Dynamic_Spectra.f90.

 --------------------------------------------------------------------
