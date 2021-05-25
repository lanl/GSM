! =============================================================================
!
!> \file
!> \brief  Contains the constructor interface for the LiteString object
!> \author CMJ (XCP-3; LANL)
!
! =============================================================================

    !> \brief Interface to the object's constructor
    public:: newLiteString
    interface newLiteString
        module procedure:: constructorMain
    end interface newLiteString

