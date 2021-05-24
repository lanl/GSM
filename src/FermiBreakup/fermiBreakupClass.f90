
! ======================================================================
!
! Module file for the simulation of fermi break up calculations.
!
! To use this physics class, users must 'use' the module, the data
! structure below (fermiFragment), the 'FermiBreakup' class, the class
! constructor (new_FermiBreakup), and 'min_fermi_A'.
!
! The parameter 'min_fermi_A' must be used to ensure that the internal
! arrays are large enough for the simulation. The arrays are recommended
! to be no smaller than 'min_fermi_A'. 'min_fermi_A" also informs users
! when Fermi Break-up physics should be used.
!
! Construction of a FermiBreakup class is described later in the
! 'new_FermiBreakup' function.
!
! To modify how much information is printed to terminal, change value of
! 'fermiVerbose' (internal only). Current available printout includes warnings,
! errors, and comments (no particle information).
!
!
! Model written and editied by many others prior to class creation.
! Written by CMJ, XCP-3, July 2018 (Initial FermiBreakup class creation).
!
! ======================================================================

module fermiBreakUpClass

! Parameters used by fermi break up routine
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use fermiBreakupParams, only : zro
  use fermiBreakupData, only : mp, mf_wnq, mpraz

  implicit none
  private

  ! For printing:
  integer(int32), private, parameter :: defaultFermiBUVerbose = 4_int32
  integer(int32), public :: fermiVerbose = defaultFermiBUVerbose

  ! Default values for options:
  integer(int32), private, parameter :: defaultRecNumNucleons = 12_int32   ! Default recommended number of nucleons by whic to simulate FBU physics
  integer(int32), private, parameter :: maxAllowedNucleons    = 16_int32   ! Max nucleus size allowed for ALL objects
  integer(int32), private, parameter :: minAllowedNucleons    =  1_int32   ! Min. nucleus size allowed for ALL objects (NOTE: nuclei w/ A=1 are simply added to the progeny bank)
  real(real64),   private, parameter :: defaultAkpScalingFlag =  0_int32   ! Flags use of some scaling factor


  ! For various results flags:
  integer(int32), public,  parameter :: noSimulationWarnings    =  0
  integer(int32), public,  parameter :: exceededProgenyArray    =  1
  integer(int32), public,  parameter :: dataApproximated        = 10
  integer(int32), public,  parameter :: fbuObjectNotConstructed = 100
  integer(int32), public,  parameter :: fbuResultNotConstructed = 200
  integer(int32), public,  parameter :: invalidNucleus          = 400


  ! For mathematical error handling:
  real(real64),   private, parameter :: div0Lim = 1.0d-15


  ! Class constructor
  public :: newFermiBreakup
  interface newFermiBreakup
     module procedure :: new_FermiBreakup
  end interface


  ! Results object constructor:
  public :: newFermiBreakUpResults
  interface newFermiBreakUpResults
     module procedure :: new_FermiBreakUpResults
  end interface


  include "fbuClassTypes.f90"


  ! Interface for the RNG and the message handling procedure:
  abstract interface
     function RANDOM() result(rndm) BIND(C)
       use, intrinsic:: iso_C_binding, only: c_double
       implicit none
       real(c_double) :: rndm
     end function RANDOM
     subroutine IOHANDLER(verbosity, type, text)
       use, intrinsic:: iso_fortran_env, only: int32
       implicit none
       integer(int32),   intent(in) :: verbosity
       integer(int32),   intent(in) :: type
       character(len=*), intent(in) :: text
     end subroutine IOHANDLER
  end interface



  ! FermiBreakup class
  type, public :: FermiBreakup
     private

     logical, private :: constructed = .FALSE.

     type(fbuIO), private :: io

     type(fermiBreakUpOptions), private :: options

     procedure(RANDOM), private, pointer, nopass :: rang => NULL()

   contains
     private
     ! Members for simulating:
     procedure, public  :: simulateFermiBU
     procedure, public  :: simulateFermiBreakUp => simulateFermiBU
     procedure, public  :: simulate             => simulateFermiBU
     procedure, public  :: disintegrate         => simulateFermiBU
     procedure, public  :: breakup              => simulateFermiBU
     procedure, public  :: deexcite             => simulateFermiBU
     procedure, public  :: decay                => simulateFermiBU
     procedure, public  :: execute              => simulateFermiBU
     procedure, public  :: start                => simulateFermiBU


     ! For querying the object and its usage
     procedure, public  :: properlyConstructed
     procedure, public  :: queryOptions
     ! (to determine appropriate FBU usage)
     procedure, private :: recommendFBUReal
     procedure, private :: recommendFBUInt
     generic,   public  :: recommendFermiBreakUp => recommendFBUInt, recommendFBUReal
     ! (Validate FBU Class Options)
     procedure, private :: validateOptions


     ! Regarding the FBU simulation:
     procedure, private :: fermid
     procedure, private :: rastar
     procedure, private :: razval
     procedure, private :: crack
     procedure, private :: divaz
     procedure, private :: wechan
     procedure, private :: tcul
     procedure, private :: disimp
     procedure, private :: wept
     procedure, private :: clpv
     procedure, private :: pint
  end type FermiBreakup


  public  :: initializeFermiData
  private :: new_FermiBreakUp
  private :: printFermi

contains

! ======================================================================

  subroutine initializeFermiData()

! ======================================================================
!
! Initialize data used by the fermi breakup simulation
!
! ======================================================================

    use fermiBreakUpData, only: initializeFermiBreakUpData
    implicit none

    call initializeFermiBreakUpData()

    return
  end subroutine initializeFermiData

! ======================================================================

  ! Allows clients to access members of the FBU object
  include "fbuAccess.f90"

  ! Constructors (and related):
  include "new_FermiBreakUp.f90"
  include "new_FermiBreakUpResults.f90"
  include "validateOptions.f90"

  ! Procedures for the simulation itself:
  include "simulateFermiBU.f90"
  include "clpv.f90"
  include "crack.f90"
  include 'disimp.f90'
  include "divaz.f90"
  include "fermid.f90"
  include "pint.f90"
  include 'rastar.f90'
  include "razval.f90"
  include "tcul.f90"
  include "wechan.f90"
  include "wept.f90"

  ! Default method for message handling:
  include "printFermi.f90"

! ======================================================================
end module fermiBreakUpClass
