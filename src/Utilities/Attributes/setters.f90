! =============================================================================
!
!> \file
!> \brief  Contains the Attribute object setter methods
!> \author CMJ (XCP-3; LANL)
!
! =============================================================================


! =============================================================================
!
!> \fn    setConstructed
!> \brief Sets the construction state of the object
!
! ARGUMENTS:
!> \param[inout] this    The Attribute (or child of) object
!> \param[in   ] flag    The new value of the construction flag
!
! =============================================================================
subroutine setConstructed(this, flag)
    implicit none
    class(Attribute), intent(inout):: this
    logical,          intent(in   ):: flag

    this%b_constructed = flag
    return
end subroutine setConstructed


! =============================================================================
!
!> \fn    setName
!> \brief Sets the name of the object
!
!> Sets the name of the Attribute (or child of) object.
!
! ARGUMENTS:
!> \param[inout] this    the attribute (or child of) object
!> \param[in   ] name    name assigned to the object
!
! =============================================================================
subroutine setName(this, name)
    implicit none
    class(Attribute), intent(inout):: this
    character(:),     intent(in   ), allocatable:: name

    this%b_name = trim(adjustl(name))
    return
end subroutine setName


! =============================================================================
!
!> \fn    setDescription
!> \brief Sets the description of the Attribute (or child of) object
!
! ARGUMENTS:
!> \param[inout] this         The attribute (or child of) object
!> \param[in   ] description  Description of the object 
!
! =============================================================================
subroutine setDescription(this, description)
    implicit none
    class(Attribute), intent(inout):: this
    character(:),     intent(in   ), allocatable:: description

    this%b_description = trim(adjustl(description))
    return
end subroutine setDescription


! =============================================================================
!
!> \fn    setLogger
!> \brief Sets the Attribute (or child of) object's Logger
!
! ARGUMENTS:
!> \param[inout] this         The attribute (or child of) object
!> \param[in   ] logger       The new Logger object
!
! =============================================================================
subroutine setLogger(this, logger)
    use Loggers, only: IOLogger => Logger
    implicit none
    class(Attribute), intent(inout):: this
    class(IOLogger),  intent(in   ), target:: logger

    this%b_log => logger
    return
end subroutine setLogger

