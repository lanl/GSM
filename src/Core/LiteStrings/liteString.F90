! =============================================================================
!>
!> \file
!> \brief  Contains the implementation of the LiteString object
!> \author CMJ
!>
! =============================================================================
!
! DATA TYPE DESCRIPTION
!
!> \class LiteString
!> \brief String object with limited functionality
!>
!> Defines a string-like object, similar to that of C++, with limited
!> functionality. The object is meant to provide some basic features that can
!> help to efficiently create strings for I/O purposes.
!
! =============================================================================
    type, public :: LiteString
        private

        !> \brief The data of the string
        character(:), private, allocatable :: body

    contains
        private

        ! >>> Alteration procedures

        !> \brief Appends text to the end of a string (concatenation)
        procedure, public :: append
        generic, public :: concat => append
        
        !> \brief Prepends text to the string, such as to add a label
        procedure, public :: prepend

        !> \brief Splits the string at the specified column
        procedure, public :: splitNearColumn

        !> \brief Erases all data in the string
        procedure, public :: clear
        generic, public :: empty => clear
        generic, public :: wipe => clear

        !> \brief Writes the output to the specified unit and clears the string
        procedure, public :: dump


        ! >>> Query procedures

        !> \brief Returns the length of the string
        procedure, public :: length

        !> \brief Returns the text contained by the string
        procedure, public :: text

    end type LiteString

