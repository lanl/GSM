
! ==============================================================================
!
! Module for the photon event generator (originally from Russia/USSR)
!
! This module contains the photon event generator class. The class:
!    1. Interpolates photon cross section data based on course angle mesh (see qintxs)
!    2. Obtains photon cross sections for "gamma + A" reactions.
!    3. Performs lorentz transformations
!
! NOTE: The object only needs to be constructed when using the interpolation!
!       The object does not use special data when obtaining photon cross section
!       or when performing the lorentz transformation.
!
!
! Written originally by others (and editied)
! Migrated to F2008 standard and a module by CMJ, XCP-3 (12/2018)
!
! ==============================================================================

module photonEventGeneratorClass

  use, intrinsic :: iso_fortran_env, only: int32, real64

  implicit none
  private

  ! For printing:
  integer(int32), private, parameter :: defaultPhotoEGVerbose = 4_int32
  integer(int32), public :: photoEGVerbose = defaultPhotoEGVerbose

  ! Constants
  real(real64), private, parameter :: zro   = 0.0_real64
  real(real64), private, parameter :: one   = 1.0_real64
  real(real64), private, parameter :: two   = 2.0_real64
  real(real64), private, parameter :: thr   = 3.0_real64
  real(real64), private, parameter :: four  = 4.0_real64
  real(real64), private, parameter :: eight = 8.0_real64
  ! (for divide by zero error protection)
  real(real64), private, parameter :: divZerLim = 1.0d-20


  public :: newPhotonEventGenerator
  interface newPhotonEventGenerator
     module procedure :: photonEGMainConstructor
  end interface


  ! Procedure pointer interface for message handling
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


  ! Provide a default I/O handling procedure for the entire module
  type, private :: photoIO
     private
     character(LEN=512),   private                  :: message =  ""
     procedure(IOHANDLER), private, pointer, nopass :: print   => printPhotoEG
  end type photoIO


  ! Data for qintxs routine (based on angle bins)
  type, private :: qintXSData
     private
     integer(int32), private :: numElements = 0_int32
     real(real64),   private, allocatable, dimension(:) :: theta
     real(real64),   private, allocatable, dimension(:) :: d
     real(real64),   private, allocatable, dimension(:) :: x11
     real(real64),   private, allocatable, dimension(:) :: x22
     real(real64),   private, allocatable, dimension(:) :: x33
     real(real64),   private, allocatable, dimension(:) :: x12
     real(real64),   private, allocatable, dimension(:) :: x23
     real(real64),   private, allocatable, dimension(:) :: x31
  end type qintXSData


  type, public :: PhotonEventGenerator
     private
     logical,          private :: constructed = .FALSE.
     type(photoIO),    private :: io
     type(qintXSData), private :: data
   contains
     private
     procedure, public  :: properlyConstructed   ! Returns the 'constructed' flag of the class
     procedure, public  :: lorentzTransform    ! Lorentz transformation of momentem to another
     procedure, public  :: photoCrossSection   ! Returns cross section based on interaction (isotopic invariance assumed)
     procedure, public  :: qintxs
     ! (called by 'photoCrossSection')
     procedure, private :: skosgdr   ! Cross section for Giant Dipole Resonance
     procedure, private :: skosga    ! Cross section for "gamma + X" reaction
     procedure, private :: skosgd    ! Cross section for "gamma + d" reaction
!     procedure, private :: skosgp    ! Cross section for  "gamma + n/p" reaction
     ! (used to establish qintxs data)
     procedure, private :: setupQintsData
  end type PhotonEventGenerator


  ! PhotonEventGenerator constructor:
  private :: photonEGMainConstructor


  private :: skosgp    ! Cross section for  "gamma + n/p" reaction
 ! (called by the skos[gP/gD/gA/GDR] procedures)
  private :: hpa
  private :: udel
  private :: wdel
  private :: wha

  ! (handles all I/O for the module)
  private :: printPhotoEG

contains

! ==============================================================================

  ! PhotonEG class constructor
  include "photonEGMainConstructor.f90"
  include "properlyConstructed.f90"

  ! Lorentz transformation
  include "lorentzTransform.f90"

  ! For obtaining gamma + A cross sections (A=1,2,etc.)
  include "photoCrossSection.f90"   ! Interface to determine the reaction
  include "skosgP.f90"                  ! gamma + P cross section
  include "skosgD.f90"                  ! gamma + D cross section
  include "skosgA.f90"                  ! gamma + A cross section
  include "skosGDR.f90"                 ! gamma + A (in GDR region) cross section
  include "skosgSupport.f90"        ! Includes those functions called by the skosgX procedures

  ! Interpolation of gamma+N differential cross sections
  include "setupQintsData.f90"
  include "qintxs.f90"

  ! Include I/O handling for the module
  include "printPhotoEG.f90"

! ==============================================================================
end module photonEventGeneratorClass
