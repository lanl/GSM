
module defaults_mod

! ==============================================================================
!
! This module file was created by CMJ, XCP-3 (8/17) with the intent to allow
! developers of GSM to more easily update GSM, to
! control what statements are printed to terminal (verbosity), and to aid in the
! management of several variables used within multiple source code files.  Such
! things include counting statistics in the CEM fashion or LAQGSM fashion,
! tallying the projectile information, using the CEM/LAQGSM RNG, etc.
!
!
! Created by CMJ, 8-25-17, XCP-3
! Modified by CMJ, 09-XX-17, XCP-3
! Modified by CMJ, 10-XX-17, XCP-3
!
! ==============================================================================

! Ensure all variables are private and declared unless otherwise noted
  use, intrinsic:: iso_fortran_env, only: int32
  implicit none
  private


! ------------------------------------------------------------------------------
! General Numbers (private)
! ------------------------------------------------------------------------------
  integer(int32), private, parameter :: zro  = 0.d0
  integer(int32), private, parameter :: one  = 1.d0
  integer(int32), private, parameter :: two  = 2.d0
  integer(int32), private, parameter :: three= 3.d0
  integer(int32), private, parameter :: four = 4.d0
  integer(int32), private, parameter :: five = 5.d0
  integer(int32), private, parameter :: six  = 6.d0
  integer(int32), private, parameter :: svn  = 7.d0
  integer(int32), private, parameter :: eight= 8.d0
  integer(int32), private, parameter :: nine = 9.d0
  integer(int32), private, parameter :: ten  = 10.d0
! ------------------------------------------------------------------------------



! ------------------------------------------------------------------------------
! For use with verbose output - when messages are printed to terminal
!      NOTE: I/O time is EXPENSIVE (computationally) - use sparingly
! =0 : NO printouts (when controlled by verbosity) except for errors
! =1 : MINIMAL print out; ONLY "NC=..."
! =2 : Prints the (1) items, and Error protection statements (when transferred)
! =3 : Prints the (2) items, and also prompts variable changes, calc. comments, and error increments
! =4 : Prints the (3) items, and also ALL WARNING/ERROR statements when re-calculating and input vars
! =5 : Prints the (4) items, and also EVERY time recalculating due to Z=0 of residuals and when parz,
!      spt exceeded in LAQhelp file
! =6 : Prints the (5) items, and EVERY error checking attempt - NOT RECOMMENDED!
!      Recommended Value: =2 (as standalone), =0 or 1 (in MCNP), =2 or 3 (when develpoing)
! ------------------------------------------------------------------------------
  integer(int32), public, parameter :: verbose = two


! ------------------------------------------------------------------------------




! ------------------------------------------------------------------------------
! Used with the Dubna IntraNuclear Cascade Model (DCM) from LAQGSM
! ------------------------------------------------------------------------------


! ----------
! Combination of cem_statistics and casc=1 (easier for GSM to determine CEM statistic use w/ DCM)
!      Recommended Value: DOES NOT MATTER - variable changed in "gsm.f90" file
  logical, public, protected :: use_DCM = .FALSE.
! ----------


! ----------
! Used to FORCE testing of the DCM.
!      Recommend Value: .FALSE.   [Only to be used when testing DCM!]
! ----------
  logical, public, parameter :: test_DCM = .FALSE.
  logical, public, parameter :: test_CEM = .FALSE.
! ----------


! ----------
! Tallies cascade particles in "upacow" (if true), in "cascaw" otherwise
!      Recommended Value: TRUE
  logical, public, parameter :: upacowTally = .TRUE.
! ----------


! ----------
! Variable controls when the spt, parz, and sptg arrays are reset
!   i.e. when there are close to 150 particle types stored, reset prior to
!        reaching this point (saftey net)
!      Recommended Value: 140
! ----------
  integer(int32), public, parameter :: resetIndx = 145;
! ----------


! ----------
! Variable is the limit of the spt, parz, and sptg arrays.  Cannot access
! anything beyond this value
!      Recommended Value: 150 (unless spt, parz, and sptg expanded)
! ----------
  integer(int32), public, parameter :: arrayLim = 150;
! ----------


! ----------
! Variable determines if unphysical/bad residuals can be used in the calculation
!      Recommended Value: FALSE
! ----------
  logical, public, parameter :: allowBadResidual = .FALSE.
! ----------


! ------------------------------------------------------------------------------


! ==============================================================================
! ------------------------------------------------------------------------------
! Declaring public subroutines
! ------------------------------------------------------------------------------
  public :: set_DCM
! ------------------------------------------------------------------------------

contains

! ----------------------------------------------------------

  subroutine set_DCM( newValue )

! ------------------------------------------------------------------------------
! This subroutine is used to set the value of use_DCM.  Routine is meant to be
! extremely simple, more of an organizational tool for the code.
! ------------------------------------------------------------------------------

    implicit none
    logical, intent(in) :: newValue

! --------------------------------------

    use_DCM = newValue

    ! Printing to user if print desired
    if ( verbose >= 3 ) then
       if ( use_DCM ) then
          write(*,1200) "DCM will be used."
       else
          write(*,1200) "DCM will NOT be used."
       endif
    endif

    return

! --------------------------------------
1200 format (3x, "Comment: ", A)
! --------------------------------------

  end subroutine set_DCM

! ==============================================================================
end module defaults_mod
