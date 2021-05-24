! =============================================================================
!
!> \file
!> \brief  Contains the constructor interface for the Attribute object
!> \author CMJ (XCP-3; LANL)
!
! =============================================================================

    !> \brief Interface to the object's constructor
    public:: newAttribute
    interface newAttribute
        module procedure:: constructorMain
    end interface newAttribute

