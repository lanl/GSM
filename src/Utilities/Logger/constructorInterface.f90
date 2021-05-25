! =============================================================================
!
!> \file
!> \brief  Contains the constructor interface for the Logger object
!> \author CMJ (XCP-3; LANL)
!
! =============================================================================

    !> \brief Interface to the object's constructor
    public:: newLogger
    interface newLogger
        module procedure:: constructorMain
    end interface newLogger

