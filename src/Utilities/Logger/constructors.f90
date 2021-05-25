! =============================================================================
!
!> \file
!> \brief  Contains the Logger object constructors
!> \author CMJ (XCP-3; LANL)
!
! =============================================================================

! =============================================================================
!
!> \fn    constructorMain
!> \brief Main Logger constructor
!
!> Main constructor for the Logger object. Clients may use this to
!> establish the level, stream, filters, etc. of the Logger object.
!
! ARGUMENTS:
!> \param[in   ] level         optional; level of message the Logger will print
!> \param[in   ] stream        optional; the stream used to print to
!> \param[in   ] filters       optional; array of message types to filter
!
! =============================================================================
 function constructorMain(level &
         & ,stream &
         & ,filters &
         & ) &
         & result(this)

     use, intrinsic:: iso_fortran_env, only: int32

     implicit none
     integer(int32), intent(in   ), optional:: level
     integer(int32), intent(in   ), optional:: stream
     integer(int32), intent(in   ), optional, dimension(:):: filters
     type(Logger):: this

! =============================================================================

     ! Interface to establish internal values of the object
     if (present(level)) call this%setLevel(level)
     if (present(stream)) call this%setStream(stream)
     if (present(filters)) call this%addFilters(filters)

     return
! =============================================================================
end function constructorMain


