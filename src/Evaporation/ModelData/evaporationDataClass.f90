
! ======================================================================
!
! This module contains the 'EvaporationData' type that is needed by
! the evaporation object. As part of the contract, users of the Evaporation
! class must first establish an 'EvaporationData' class, which is then
! passed in to the Evaporation object's constructor.
!
! This data is NOT to be modified by the user, however must be accessible
! to the evaporation class.
!
! Written by CMJ, XCP-3, 8/2018 (Original Class Creation)
!
! ======================================================================

module evaporationDataClass

  use, intrinsic:: iso_fortran_env, only: int32, real64
  use evaporationFissionData, only: iiz, inn, levelDensityGCCIFlag

  implicit none
  private


  ! For verbosity:
  integer(int32), private, parameter :: defaultEvapDataVerbose = 4_int32
  integer(int32), public :: evapDataVerbose = defaultEvapDataVerbose

  ! Level density parameter option; =0 (GCCI), =1 (a=A/8), >1 (a=a/alev)
  real(real64),   public,  parameter :: defaultAlev = levelDensityGCCIFlag


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


  ! To construct evaporation data class
  public :: newEvaporationData
  interface newEvaporationData
     module procedure :: evapDCMainConstructor
  end interface

  type, private :: evaporationDataVariables
     private
     real(real64), private, dimension(iiz, inn) :: fact10 = 0_real64
     real(real64), private, dimension(iiz, inn) :: sx0    = 0_real64
     real(real64), private, dimension(iiz, inn) :: taux0  = 0_real64
     real(real64), private, dimension(iiz, inn) :: ax0    = 0_real64
     real(real64), private, dimension(iiz, inn) :: tx0    = 0_real64
     real(real64), private, dimension(iiz, inn) :: extx0  = 0_real64
  end type evaporationDataVariables


  type, private :: evaporationDataIO
     private
     character(LEN=512),   private                  :: message =  ""
     procedure(IOHANDLER), private, pointer, nopass :: print   => printEvaporation
  end type evaporationDataIO


  type, public  :: evaporationDataOptions
     private
     real(real64), public :: alev = defaultAlev
  end type evaporationDataOptions


  ! Evaporation data class
  type, public :: EvaporationData
     private

     logical, private :: constructed = .FALSE.
     type(evaporationDataIO),        private :: io
     type(evaporationDataOptions),   private :: options
     type(evaporationDataVariables), private :: data

   contains
     private
     ! To verify construction:
     procedure, public  :: properlyConstructed

     ! Access to data members:
     procedure, public  :: alev
     procedure, public  :: fact10
     procedure, public  :: sx0
     procedure, public  :: taux0
     procedure, public  :: ax0
     procedure, public  :: tx0
     procedure, public  :: extx0
     procedure, private :: checkIndex
     procedure, public  :: queryOptions
     procedure, private :: validateOptions
  end type EvaporationData


contains

! ======================================================================

  ! Class constructor
  include "evapDCMainConstructor.f90"
  include "validateOptions.f90"

  include "geta.f90"


  ! Includes client access functions (only "getters"):
  include "evapDataAccess.f90"


  ! For printing messages
  include "printEvaporation.f90"


  function properlyConstructed ( data ) result(constructed)

! ======================================================================
!
! This procedure returns to the client a logical flag indicating if the
! data object was properly constructed or not
!
! ======================================================================

    implicit none
    class(EvaporationData), intent(in   ) :: data
    logical                               :: constructed

! ======================================================================

    constructed = data%constructed

    return
! ======================================================================
  end function properlyConstructed


! ======================================================================
end module evaporationDataClass
