
! ==============================================================================
!
! This file contains information for the number of times each model was used
!
!
! Written by CMJ, XCP-3, 01/2019
!
! ==============================================================================

  type, private :: gsmModelUsage
     private
     ! Number of times a model was simulated
     integer(int64), public  :: numSDCMSim     = 0_int64
     integer(int64), public  :: numMDCMSim     = 0_int64
     integer(int64), public  :: numCoalesSim   = 0_int64
     integer(int64), public  :: numPreeqSim    = 0_int64
     integer(int64), public  :: numEvapSim     = 0_int64
     integer(int64), public  :: numFermiBUSim  = 0_int64


     ! For sDCM/mDCM specifically:
     integer(int64), public  :: numSDCMRestarts = 0_int32
     integer(int64), public  :: numMDCMRestarts = 0_int32


     ! Number of events ending from a specific physics (Evaporation most likely)
     !    NOTE: This could also be due to the progeny array being exceeded
     !    NOTE: An event should NEVER end in coalescence (it reduces the number of progeny)
     integer(int64), public  :: numSDCMEnd     = 0_int64
     integer(int64), public  :: numMDCMEnd     = 0_int64
     integer(int64), public  :: numCoalesEnd   = 0_int64
     integer(int64), public  :: numPreeqEnd    = 0_int64
     integer(int64), public  :: numEvapEnd     = 0_int64
     integer(int64), public  :: numFermiBUEnd  = 0_int64
  end type gsmModelUsage
