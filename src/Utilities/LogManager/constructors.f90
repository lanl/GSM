! =============================================================================
!
!> \file
!> \brief  Contains the LogManager object constructors
!> \author CMJ (XCP-3; LANL)
!
! =============================================================================

! =============================================================================
!
!> \fn    constructorMain
!> \brief Main LogManager constructor
!
!> Main constructor for the LogManager object. Clients may use this to
!> establish the add Logger objects to print messages to certain streams
!
! ARGUMENTS:
!> \param[in   ] loggerArray    optional; array of Logger objects to provide
!>
!
! =============================================================================
 function constructorMain(loggerArray &
         & ) &
         & result(this)

     use, intrinsic:: iso_fortran_env, only: int32
     use Loggers, only: Logger

     implicit none
     class(Logger),  intent(in   ), optional, dimension(:):: loggerArray
     type(LogManager):: this

     integer(int32):: indx

! =============================================================================

     ! Interface to establish internal values of the object
     if (present(loggerArray)) then
         if (size(loggerArray) > 0) then
             call this%resize (size(loggerArray))
             addLoop: do indx = 1_int32, size(loggerArray)
                 call this%addLogger (loggerArray(indx))
             end do addLoop
         end if
     end if

     return
! =============================================================================
end function constructorMain


