! =============================================================================
!
!> \file
!> \brief  Contains the abstraction for OpenACC
!> \author CMJ (XCP-3; LANL)
!>
!> The module provides the abstraction layer for calls to the OpenACC library
!> when it is defined (OpenACC) by the compiler.
!>
!> Note that OpenACC funcionality is only invoked when the project is compiled
!> with the \c -fopenacc compiler flag via the \c _OPENACC preprocessor definition.
!
! =============================================================================
module OpenACCAbstraction

#if defined(_OPENACC)
  use OPENACC
#endif

  ! Import Contracts
  use Contracts

  implicit none
  private

  ! Public functions
  public :: &
       & gpuThreadingAvailable, &
       & numCurrentDevices

contains


! =========================================================================== !
!
!> \fn
!> \brief Returns a logical flag indicating if OpenACC is available
!
! =========================================================================== !
  function gpuThreadingAvailable() result(available)
    implicit none
    logical :: available
#if defined(_OPENACC)
    available = .TRUE.
#else
    available = .FALSE.
#endif
    return
  end function gpuThreadingAvailable


! =========================================================================== !
!
!> \fn
!> \brief Returns the number of devices currently in action

! =========================================================================== !
  function numCurrentDevices() result(numDevices)
    use, intrinsic:: iso_fortran_env, only: int32
    implicit none
    integer(int32) :: numDevices
#if defined(_OPENACC)
    numDevices = acc_get_num_devices(acc_get_device_type())
#else
    numDevices = 1_int32
#endif
    return
  end function numCurrentDevices


! =========================================================================== !
end module OpenACCAbstraction
