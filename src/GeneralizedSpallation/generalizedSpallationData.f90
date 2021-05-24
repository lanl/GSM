
! ==============================================================================
!
! This module contains the data necessary for the GSM class:
!
! ==============================================================================

module generalizedSpallationData

  use, intrinsic:: iso_fortran_env, only: int32, int64, real64
  use randomNumberGenerator, only: rang

  implicit none
  private

  ! Flags if data was constructed:
  logical, public, protected :: gsmDataInitialized

  ! For message filtering:
  integer(int32), private, parameter :: defaultGSMVerbose = 4_int32
  integer(int32), public :: gsmVerbose = defaultGSMVerbose

  ! For the random number generator:
  ! NOTE: A (0) value defers the value to its default within the module for the specified type
  integer(int32), private, parameter :: defaultRndmType       = 1_int32
  integer(int64), private, parameter :: defaultSeed           = 0_int64
  integer(int64), private, parameter :: defaultStride         = 0_int64
  integer(int64), private, parameter :: defaultNumSeedAdvance = 0_int64
  integer(int32), private, parameter :: defaultPrintInfo      = 1_int32   ! Flags that info will be printed


  ! Message handling procedure and random number generator:
  abstract interface
     subroutine IOHANDLER(verbosity, type, text)
       use, intrinsic:: iso_fortran_env, only: int32
       implicit none
       integer(int32),   intent(in) :: verbosity
       integer(int32),   intent(in) :: type
       character(len=*), intent(in) :: text
     end subroutine IOHANDLER
     function RANDOM() result(rndm)
       use, intrinsic:: iso_fortran_env, only: real64
       implicit none
       real(real64) :: rndm
     end function RANDOM
  end interface


  ! Message handler:
  type, private :: gsmDataIO
     private
     character(LEN=512), private :: message = ""
     procedure(IOHANDLER), private, nopass, pointer :: print => printGSM
  end type gsmDataIO
  type(gsmDataIO), private :: dataIO


  ! Default pointer for the RNG:
  procedure(RANDOM), private, pointer :: rndmProcedure => rang


  include "gsmData1.f90"


  public  :: initializeGSMData
  public  :: queryRNG
  private :: initializeSubModelData
  private :: printGSM

contains

! ==============================================================================

  function queryRNG() result(rng)
    ! Return the random number generator pointer to the calling routine
    implicit none
    procedure(RANDOM), pointer :: rng
    rng => rndmProcedure
    return
  end function queryRNG

! ==============================================================================

  include "initializeGSMData.f90"

  include "initializeSubModelData.f90"

  include "printGSM.f90"

! ==============================================================================
end module generalizedSpallationData
