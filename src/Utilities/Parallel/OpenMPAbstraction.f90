! =============================================================================
!
!> \file
!> \brief  Contains the abstraction for OpenMP
!> \author CMJ (XCP-3; LANL)
!>
!> The module provides the abstraction layer for calls to the OpenMP library
!> when it is defined (OpenMP) by the compiler. \b Note: This abstraction layer
!> is based on that used in the Nemesis C++ environment hosted by the Shift
!> team (within SCALE) at ORNL. Note that OpenMP funcionality is only invoked
!> when the project is compiled with the \c -fopenmp compiler flag via the
!> \c _OPENMP preprocessor definition.
!
! =============================================================================
module OpenMPAbstraction

#if defined(_OPENMP)
  use OMP_LIB
#endif

  ! Import Contracts
  use Contracts

  implicit none
  private

  ! Public functions
  public :: &
       & multiThreadingAvailable, &
       & cancellationEnabled, &
       & setNumThreads, &
       & numCurrentThreads, &
       & threadID, &
       & numAvailableThreads, &
       & numAvailableProcessors, &
       & useDynamicThreading, &
       & dynamicThreading, &
       & threadTime

contains


! =========================================================================== !
!
!> \fn
!> \brief Returns a logical flag indicating if OpenMP is available
!
! =========================================================================== !
  function multiThreadingAvailable() result(available)
    implicit none
    logical :: available
#if defined(_OPENMP)
    available = .TRUE.
#else
    available = .FALSE.
#endif
    return
  end function multiThreadingAvailable


! =========================================================================== !
!
!> \fn
!> \brief Query if OpenMP cancellation policy is enabled

! =========================================================================== !
  function cancellationEnabled() result(enabled)
    implicit none
    logical :: enabled
#if defined(_OPENMP)
    enabled = omp_get_cancellation()
#else
    enabled = .FALSE.
#endif
    return
  end function cancellationEnabled


! =========================================================================== !
!
!> \fn
!> \brief Sets the number of threads to use in the current execution

! =========================================================================== !
  subroutine setNumThreads(numThreads)
    use, intrinsic:: iso_fortran_env, only: int32
    implicit none
    integer(int32), intent(in   ) :: numThreads
#if defined(_OPENMP)
    call insist(numThreads > 0, &
        & "Number of threads cannot be less than 0.", &
        & __FILE__, __LINE__)
    call omp_set_num_threads(numThreads)
#endif
    return
  end subroutine setNumThreads


! =========================================================================== !
!
!> \fn
!> \brief Returns the number of threads currently in action

! =========================================================================== !
  function numCurrentThreads() result(numThreads)
    use, intrinsic:: iso_fortran_env, only: int32
    implicit none
    integer(int32) :: numThreads
#if defined(_OPENMP)
    numThreads = omp_get_num_threads()
#else
    numThreads = 1_int32
#endif
    return
  end function numCurrentThreads


! =========================================================================== !
!
!> \fn
!> \brief Returns the thread ID

! =========================================================================== !
  function threadID() result(ID)
    use, intrinsic:: iso_fortran_env, only: int32
    implicit none
    integer(int32) :: ID
#if defined(_OPENMP)
    ID = omp_get_thread_num()
#else
    ID = 1_int32
#endif
    return
  end function threadID


! =========================================================================== !
!
!> \fn
!> \brief Returns the number of threads available in a thread block
!>
!> Returns the number of threads available in a thread block. This \b always
!> returns the number of threads set using \c set_num_threads().

! =========================================================================== !
  function numAvailableThreads() result(numAvail)
    use, intrinsic:: iso_fortran_env, only: int32
    implicit none
    integer(int32) :: numAvail
#if defined(_OPENMP)
    numAvail = omp_get_max_threads()
#else
    numAvail = 1_int32
#endif
    return
  end function numAvailableThreads


! =========================================================================== !
!
!> \fn
!> \brief Returns the number of processors on the device

! =========================================================================== !
  function numAvailableProcessors() result(numAvailable)
    use, intrinsic:: iso_fortran_env, only: int32
    implicit none
    integer(int32) :: numAvailable
#if defined(_OPENMP)
    numAvailable = omp_get_num_procs()
#else
    numAvailable = 1_int32
#endif
    return
  end function numAvailableProcessors


! =========================================================================== !
!
!> \fn
!> \brief Turns on or off dynamic adjustment of the number of threads

! =========================================================================== !
  subroutine useDynamicThreading(isDynamic)
    implicit none
    logical, intent(in   ) :: isDynamic
#if defined(_OPENMP)
    call omp_set_dynamic(isDynamic)
#endif
    return
  end subroutine useDynamicThreading


! =========================================================================== !
!
!> \fn
!> \brief Returns if dynamic adjustment of the available threads is allowed

! =========================================================================== !
  function dynamicThreading() result(dynamicThread)
    implicit none
    logical :: dynamicThread
#if defined(_OPENMP)
    dynamicThread = omp_get_dynamic()
#else
    dynamicThread = .FALSE.
#endif
    return
  end function dynamicThreading


! =========================================================================== !
!
!> \fn
!> \brief Returns the elapsed wall clock time in seconds for a thread

! =========================================================================== !
  function threadTime() result(time)
    use, intrinsic:: iso_fortran_env, only: real64
    implicit none
    real(real64) :: time
#if defined(_OPENMP)
    time = omp_get_wtime()
#else
    time = 0.0_real64
#endif
    return
  end function threadTime


! =========================================================================== !
end module OpenMPAbstraction
