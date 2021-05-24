
! ====================================================================
!
! Module file containing the class that simulated the physics of the 
! Standard Dubna Cascade Model.
!
!
! Written by CMJ, XCP-3, 12/2018
!
! ====================================================================
module standardDCMClass

  use, intrinsic:: iso_fortran_env, only: int32, real64
  use standardDCMParams, only: zro, one
  use standardDCMDataClass, only: standardDCMData
  use molnixClass, only: Molnix


  implicit none
  private

  ! For printing:
  integer(int32), private, parameter :: defaultSDCMVerbose = 4_int32
  integer(int32), public :: sDCMVerbose = defaultSDCMVerbose

  ! Default behaviors
  ! (regarding inclusion of refraction and reflection at nuclear boundaries; =0 if OFF)
  integer(int32), private, parameter :: boundaryEffectsDefault = 0_int32
  ! (regarding use of new angular distribution approximations) - new ones correct the unphysical "pole-like" behavior
  integer(int32), private, parameter :: newAngularDistributionsDefault = 1_int32


  ! For divide by zero checks
  real(real64), private, parameter :: div0Lim = 1.0d-15


  ! Interface to establish a new S-DCM class
  public :: newStandardDCM
  interface newStandardDCM
     module procedure :: sDCMMainConstructor
  end interface


  ! Interface to establish a new results objecty for the S-DCM class
  public :: newStandardDCMResults
  interface newStandardDCMResults
     module procedure :: new_StandardDCMResults
  end interface


  ! Procedure pointer interface
  abstract interface
     ! Interface for I/O handling
     subroutine IOHANDLER(verbosity, type, text)
       use, intrinsic:: iso_fortran_env, only: int32
       implicit none
       integer(int32),   intent(in) :: verbosity
       integer(int32),   intent(in) :: type
       character(len=*), intent(in) :: text
     end subroutine IOHANDLER
     function RANDOM() result(rndmNum) BIND(C)
       use, intrinsic:: iso_C_binding, only: c_double
       implicit none
       real(c_double) :: rndmNum
     end function RANDOM
  end interface

  ! Include target and projectile data type
  include "sDCMClassTypes.f90"

  type, public :: StandardDCM
     private

     ! Flag for proper construction
     logical, private :: constructed = .FALSE.

     ! I/O Handling
     type(sDCMIO), private :: io

     ! Options for the simulation
     type(sDCMOptions),  private :: options

     ! Data for the sDCM Class
     type(Molnix),          private, pointer :: molnixE  => NULL()

     ! Random number generator
     procedure(RANDOM), private, nopass, pointer :: rang => NULL()

     ! NOTE: Results and calculation data types are created in the main routine,
     !       minimizing the state of the class
   contains
     private
     ! Public  Member(s)
     !    (Various procedure names to use to simulate)
     procedure, public :: simulate         => cascad
     procedure, public :: execute          => cascad
     procedure, public :: start            => cascad
     procedure, public :: interact         => cascad
     procedure, public :: collide          => cascad
     procedure, public :: simulateCascade  => cascad
     procedure, public :: simulateINC      => cascad
     !    (Accessing members of the class)
     procedure, public :: properlyConstructed
     procedure, public :: queryOptions
     procedure, public :: queryMolnix
     procedure, public :: queryRNG
     ! ------------------
     ! Private Member(s)
     procedure, private :: abel
     procedure, private :: absorp
     procedure, private :: bindnuc
     procedure, private :: binel
     procedure, private :: cascad
     procedure, private :: cduarte
     procedure, private :: chabs
     procedure, private :: chinel
     procedure, private :: cosel
     procedure, private :: cosgamn
     procedure, private :: costa
     procedure, private :: costan
     procedure, private :: cosex
     procedure, private :: direct8
     procedure, private :: elex
     procedure, private :: geom8
     procedure, private :: isobar
!     procedure, private :: jtypa
!     procedure, private :: jtypb
     procedure, private :: kinema
     procedure, private :: paulip
     procedure, private :: partn
     procedure, private :: pinpn
     procedure, private :: pinpn1
!     procedure, private :: pmom
     procedure, private :: pointe
     procedure, private :: pointe1
!     procedure, private :: poten
     procedure, private :: qints
     procedure, private :: refrac
     procedure, private :: refrac1
     procedure, private :: rotor
     procedure, private :: setPhotoChannelSigma
     procedure, private :: sigmat8
!     procedure, private :: slqek
     procedure, private :: statSDCM
     procedure, private :: tinvu
     procedure, private :: typint
     procedure, private :: vmnsp
     procedure, private :: wim
     procedure, private :: wopt
  end type StandardDCM



  ! Non-Standard DCM routines (not in class)
  private :: sDCMMainConstructor         ! Gives client a SDCM data class
  private :: new_StandardDCMResults      ! Gives client a "results" object to pass in to sDCM Class
  private :: printSDCM                   ! Handles all printing from SDCM class (default)

  ! Used by sDCM simulation (but not needed to be part of the class):
  private :: slqek
  private :: poten
  private :: pmom
  private :: jtypa
  private :: jtypb

contains

! ====================================================================

  ! Allow client to check if sDCM was properly constructed
  include "sDCMAccess.f90"


  ! Include sDCM constructor functions
  include "sDCMMainConstructor.f90"


  ! Returns a standard DCM results type
  include "new_StandardDCMResults.f90"



  ! Default I/O Handler for the data class
  include "printSDCM.f90"



  ! Import all other procedures (internal to SDCM Object)
  include "abel.f90"
  include "absorp.f90"
  include "bindnuc.f90"
  include "binel.f90"
  include "cascad.f90"
  include "cduarte.f90"
  include "chabs.f90"
  include "chinel.f90"
  include "cosel.f90"
  include "costa.f90"
  include "direct.f90"
  include "elex.f90"
  include "geom.f90"
  include "isobar.f90"
  include "jtypa.f90"
  include "jtypb.f90"
  include "kinema.f90"
  include "partn.f90"
  include "paulip.f90"
  include "pinpn.f90"
  include "pmom.f90"
  include "pointe.f90"
  include "poten.f90"
  include "qints.f90"
  include "refrac.f90"
  include "rotor.f90"
  include "setPhotoChannelSigma.f90"
  include "sigmat.f90"
  include "slqek.f90"
  include "stat.f90"
  include "tinvu.f90"
  include "typint.f90"
  include "vmnsp.f90"
  include "wim.f90"
  include "wopt.f90"

! ====================================================================
end module standardDCMClass
