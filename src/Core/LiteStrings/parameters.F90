! =============================================================================
!>
!> \file
!> \brief  Contains the LiteStrings parameters and default values
!> \author CMJ
!>
! =============================================================================
   
    !> \brief Flags to use scientific form when converting to strings
    logical, private, parameter :: useScientificForm = .False.

    !> \brief Format specifiers for different numerical conversions
    !> @{
    character(*), public, parameter :: intFormat = "(i0)"
    character(*), public, parameter :: realFormat = "(f0.3)"
    character(*), public, parameter :: scientificFormat = "(es13.4)"
    character(*), public, parameter :: logicalFormat = "(l)"
    character(*), public, parameter :: stringFormat = "(A)"
    !> @}

    !> \brief Default value at which to split a string when dumping
    integer(gsmInt8), private, parameter :: defaultColumnSplit = 85

    !> \brief Range to facilitate finding good breakpoint to split a string by
    integer(gsmInt8), private, parameter :: columnSplitRange = 10

