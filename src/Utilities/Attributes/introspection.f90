! =============================================================================
!
!> \file
!> \brief  Contains the Attribute introspection (state getters) procedures
!> \author CMJ (XCP-3; LANL)
!
! =============================================================================


! =============================================================================
!
!> \fn    getConstructed
!> \brief Returns the flagged construction state of the object
!
! ARGUMENTS:
!> \param[inout] this    Attribute (or child of) object
!
! =============================================================================
function getConstructed(this) result(constructed)
      implicit none
      class(Attribute), intent(in   ):: this
      logical:: constructed

      constructed = this%b_constructed
      return
end function getConstructed


! =============================================================================
!
!> \fn    getName
!> \brief Returns the name of the object
!
! ARGUMENTS:
!> \param[inout] this    Attribute (or child of) object
!
! =============================================================================
function getName(this) result(name)
    implicit none
    class(Attribute), intent(in   ):: this
    character(:), allocatable:: name

    name = this%b_name
    return
end function getName


! =============================================================================
!
!> \fn    getDescription
!> \brief Returns the description of the object
!
! ARGUMENTS:
!> \param[inout] this    Attribute (or child of) object
!
! =============================================================================
function getDescription(this) result(description)
    implicit none
    class(Attribute), intent(in   ):: this
    character(:), allocatable:: description

    description = this%b_name
    return
end function getDescription


! =============================================================================
!
!> \fn    getLogger
!> \brief Returns the Logger pointer of the object
!
! ARGUMENTS:
!> \param[inout] this    Attribute (or child of) object
!
! =============================================================================
function getLogger(this) result(logger)
    use Loggers, only: IOLogger => Logger
    implicit none
    class(Attribute), intent(in   ):: this
    type(IOLogger), pointer:: logger

    logger => this%b_log
    return
end function getLogger

