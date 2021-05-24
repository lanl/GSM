
! ====================================================================
!
! Module file that contains the parameterized data utilized by the
! standard Dubna Cascade Model (DCM)
!
! Written by CMJ, XCP-3, 12/2018
!
! ====================================================================

module standardDCMDataClass

  use, intrinsic:: iso_fortran_env, only: int32, real64
  use standardDCMDataParams, only: zro, one

  implicit none
  private

  ! For printing:
  integer(int32), private, parameter :: defaultSDCMDataVerbose = 4_int32
  integer(int32), public :: sDCMDataVerbose = defaultSDCMDataVerbose

  ! For construction errors:
  integer(int32), public, parameter :: objectNotConstructed = -1_int32
  integer(int32), public, parameter :: properlyConstructed  =  0_int32
  integer(int32), public, parameter :: invalidTargetFlag    =  1_int32
  integer(int32), public, parameter :: targetInitFailed     = 10_int32


  ! Default options:
  ! (default number of nuclear zones)
  integer(int32), private, parameter :: numZonesDefault =  7_int32
  integer(int32), private, parameter :: minZonesAllowed =  1_int32
  integer(int32), private, parameter :: maxZonesAllowed = 10_int32
  ! (regarding nucleon densities)
  ! rho(r) = rho0 { 1 - exp[ (r - r0 A**(1/3) ) / expDenom ] }
  real(real64),   private, parameter :: expDenomDefault =  0.545_real64   ! Value 
  real(real64),   private, parameter :: r0Default       =  1.07_real64    ! 
  real(real64),   private, parameter :: maxRadDefault   = 10.0_real64     ! 
  real(real64),   private, parameter, dimension(maxZonesAllowed) :: aveDenDefault = &
       & [  0.95_real64, 0.80_real64, 0.50_real64, 0.20_real64,  &
       &    0.10_real64, 0.05_real64, 0.01_real64, 0.00_real64,  &
       &    0.00_real64, 0.00_real64                             ]
  ! (regarding pion potentials and separation/binding energies)
  real(real64),   private, parameter :: pionPoteDefault  = 0.025_real64   ! Potential well used for pions [GeV]
  real(real64),   private, parameter :: sepEnergyDefault = 0.007_real64   ! Single nucleon binding energy [GeV]


  ! For divide by zero checks
  real(real64),  private, parameter :: div0Lim = 1.0d-15


  ! Interface to establish a new data class
  public :: newStandardDCMData
  interface newStandardDCMData
     module procedure :: sDCMDataMainConstructor
  end interface


  abstract interface
     ! Interface for I/O handling
     subroutine IOHANDLER(verbosity, type, text)
       use, intrinsic:: iso_fortran_env, only: int32
       implicit none
       integer(int32),   intent(in) :: verbosity
       integer(int32),   intent(in) :: type
       character(len=*), intent(in) :: text
     end subroutine IOHANDLER
  end interface

  ! Include target data type
  include "sDCMDataClassTypes.f90"

  type, public :: StandardDCMData
     private
     ! Flag indicating if the data object was constructed
     logical, private :: objectConstructed = .FALSE.
     integer(int32), private :: constructed = objectNotConstructed

     ! I/O Handling (all messages filter through this type)
     type(sDCMDataIO),      private :: io

     ! Contains all options for the object
     type(sDCMDataOptions), private :: options

     ! Data for the sDCM Class
     type(sDCMDataTarget),  private :: target

   contains
     private
     ! For querying object:
     procedure, public :: properlyConstructedObject
     procedure, public :: properlyConstructed => properlyConstructedObject
     procedure, public :: constructionState
     procedure, public :: queryOptions

     ! For accessing target data:
     procedure, public :: zoneBoundR
     procedure, public :: zoneBRFrac
     procedure, public :: protonDensity
     procedure, public :: neutronDensity
     procedure, public :: coulombPote
     procedure, public :: protFermiMom
     procedure, public :: neutFermiMom
     procedure, public :: numBaryons
     procedure, public :: numProtons
     procedure, public :: numZones
     procedure, public :: aTargThrd
     procedure, public :: geomCrossSection
     procedure, public :: pionPote
     procedure, public :: setSepEnergy
     procedure, public :: getSepEnergy

     ! Private Member(s)
     procedure, private :: checkIndex   ! Checks for valid element from the target array (if not, returns closest value)
     procedure, private :: setupTarget
     procedure, private :: validateOptions

  end type StandardDCMData

  abstract interface
     subroutine IntegralInterface(sxfi, bs, rs, fis, dataObj)
       use, intrinsic:: iso_fortran_env, only: real64
       import StandardDCMData
       implicit none
       real(real64), intent(in   ) :: sxfi
       real(real64), intent(in   ) :: bs
       real(real64), intent(in   ) :: rs
       real(real64), intent(  out) :: fis
       type(StandardDCMData), intent(inout) :: dataObj
     end subroutine IntegralInterface
  end interface

  ! Non-class specific routines
  private :: sDCMDataMainConstructor   ! Gives client SDCM data class
  private :: printSDCMData         ! Handles all printing from SDCM data class

  private :: sfint
  private :: sfint1
  private :: sfis
  private :: ss1
  private :: ss2

contains

! ====================================================================

  ! Include constructor and its associated procedures:
  include "sDCMDataMainConstructor.f90"
  include "validateOptions.f90"
  include "setupTarget.f90"

  ! Include public procedures for clients to access data
  include "sDCMDataAccess.f90"      ! Gives a client access to the members of'target'
  include "constructionState.f90"   ! Tells client if construction was successful

  ! Include module-specific routines (used by data class, but class use not necessary
  include "sfint.f90"
  include "sfis.f90"
  include "ss12.f90"

  ! Default I/O Handler for the data class
  include "printSDCMData.f90"

! ====================================================================
end module StandardDCMDataClass
