
! ======================================================================
!
! Module for the Molnix class.
!
! Clients may use the Molnix class by 'use'-ing the module and "Molnix"
! (i.e. "use molnixClass, only: Molnix"). An instance of the Molnix class
! may be created according to the specification listed in the constructor,
! "new_Molnix".
!
!
! Written by CMJ, XCP-3, August 2018 (Molnix class creation)
!
! ======================================================================

module molnixClass

  use, intrinsic :: iso_fortran_env, only: int32, real64

  implicit none
  private
  ! For printing:
  integer(int32), private, parameter :: defaultMolnixVerbose = 4_int32
  integer(int32), public :: molnixVerbose = defaultMolnixVerbose


  ! Module defaults
  real(real64),   private, parameter :: defaultCevap  = 12.0_real64


  ! For divide by zero errors
  real(real64),   private, parameter :: div0Lim = 1.0d-15

  ! Load experimental data
  include "molnixData.f90"


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


  ! Class interface
  public :: newMolnix
  interface newMolnix
     module procedure :: new_Molnix
  end interface newMolnix


  ! Options type:
  type, public :: molnixOptions
     private
     real(real64), public :: cevap = defaultCevap
  end type molnixOptions


  ! I/O type:
  type, private :: molnixIO
     private
     character(LEN=512),   private                  :: message  = ""
     procedure(IOHANDLER), private, pointer, nopass :: print    => printMolnix
  end type molnixIO


  ! Molnix class
  type, public :: Molnix
     private

     logical,             private :: constructed = .FALSE.
     type(molnixOptions), private :: options
     type(molnixIO),      private :: io

   contains
     private
     ! For obtaining results
     procedure, public  :: massExcess
     procedure, public  :: shellEnergy
     procedure, public  :: pairingGap
     procedure, public  :: defineEnergy
     procedure, private :: mnmacro

     ! For obtaining the behaviors/numbers used by the class
     procedure, public  :: properlyConstructed
     procedure, public  :: queryOptions

     ! For accessing data relating to the model:
     procedure, public :: nmina
     procedure, public :: nmaxa
     procedure, public :: nmin
     procedure, public :: nmax

  end type Molnix


  private :: new_Molnix
  private :: printMolnix

contains


! ======================================================================

  ! Class constructor
  include "new_Molnix.f90"

  ! Routine that users can call
  include "shellEnergy.f90"
  include "massExcess.f90"
  include "pairingGap.f90"
  include "defineEnergy.f90"
  include "mnmacro.f90"

  ! For accessing model data and private class members
  include "dataAccess.f90"

  ! For message handling:
  include "printMolnix.f90"

! ======================================================================
end module molnixClass
