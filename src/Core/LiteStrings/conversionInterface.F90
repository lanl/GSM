! =============================================================================
!
!> \file
!> \brief  Contains the numerical conversion interface
!> \author CMJ
!
! =============================================================================

    !> \brief Conversion interface to convert numbers to a string format.
    public:: toString
    interface toString
        module procedure:: int8ToString
        module procedure:: int16ToString
        module procedure:: int32ToString
        module procedure:: int64ToString
        module procedure:: floatToString
        module procedure:: doubleToString
        module procedure:: longDoubleToString
        module procedure:: logicalToString
        module procedure:: stringToString
    end interface toString

