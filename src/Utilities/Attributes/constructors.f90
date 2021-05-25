! =============================================================================
!
!> \file
!> \brief  Contains the Attribute object constructors
!> \author CMJ (XCP-3; LANL)
!
! =============================================================================

! =============================================================================
!
!> \fn    constructorMain
!> \brief Main Attribute constructor
!
!> Main constructor for the Attribute object. Clients may use this to
!> set the name, description, Logger, etc. of the object.
!
! ARGUMENTS:
!> \param[in   ] name           name of the object
!> \param[in   ] description    optional; description of the object
!> \param[in   ] logger         optional; Logger for the object
!
! =============================================================================
 function constructorMain(name &
         & ,description &
         & ,logger &
         & ) &
         & result(this)

     use, intrinsic:: iso_fortran_env, only: int32
     use Loggers, only: MsgLogger => Logger

     implicit none
     character(:),  intent(in   ), allocatable:: name
     character(:),  intent(in   ), allocatable, optional:: description
     class(MsgLogger), intent(in   ), target,   optional:: logger
     type(Attribute):: this

     integer(int32):: allocateStatus

! =============================================================================

     ! Flag as constructed
     call this%setConstructed(.TRUE.)

     ! Interface to establish internal values of the object
     call this%setName(name)
     if (present(description)) call this%setDescription(description)
     if (present(logger)) then
         call this%setLogger(logger)
     else
         allocate(this%b_log, stat=allocateStatus)
         if (allocateStatus /= 0_int32) then
             call this%setConstructed(.FALSE.)
             !> \todo Add exception call here (NotImplemented)
         end if
     end if

     return
! =============================================================================
end function constructorMain


