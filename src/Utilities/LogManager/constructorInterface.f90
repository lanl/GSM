! =============================================================================
!
!> \file
!> \brief  Contains the constructor interface for the LogManager object
!> \author CMJ (XCP-3; LANL)
!
! =============================================================================

    !> \brief Interface to the object's constructor
    public:: newLogManager
    interface newLogManager
        module procedure:: constructorMain
    end interface newLogManager

