! =============================================================================
!
!> \file
!> \brief  Contains the implementation of the Logger object
!> \author CMJ (XCP-3; LANL)
!
! =============================================================================
!
! DATA TYPE DESCRIPTION
!
!> \class Logger
!> \brief Defines the base Logger object
!
!> Defines the base Logger object used by GSM and its sub-models and
!> model-data.
!
! =============================================================================
type, public:: Logger
    private

    !> \brief Sets the verbosity of the Logger (i.e. filters messages)
    integer(int32), private:: b_level = defaultLevel

    !> \brief The stream the Logger prints to
    integer(int32), private:: b_stream  = defaultStream

    !> \brief Filters (i.e. removes) messages of the type from the stream
    integer(int32), private, dimension(:), allocatable:: b_filters

    !> \brief Quick reference for the number of filters present
    integer(int32), private:: b_numFilters = 0_int32

 contains
    private

    ! >>> LOGGING PROCEDURES

    !> \brief Prints critical errors
    procedure, public:: critical
    generic, public:: fatal => critical

    !> \brief Prints generic errors
    procedure, public:: error

    !> \brief Prints warnings
    procedure, public:: warning
    generic, public:: warn => warning

    !> \brief Prints general info
    procedure, public:: info
    generic, public:: log => info

    !> \brief Prints debugging information
    procedure, public:: debug

    !> \brief Prints all information
    procedure, public:: print

    !> \brief Prints a message to the Logger's stream
    procedure, private:: baseLog


    ! >>> INTROSPECTION

    !> \brief Returns the message filter used
    procedure, public:: getLevel

    !> \brief Returns the stream printed to
    procedure, public:: getStream

    !> \brief Returns the number of filters
    procedure, public:: getNumFilters

    !> \brief Returns the filtered message types
    procedure, public:: getFilters

    !> \brief Returns if a message type is to be printed or not
    procedure, private:: printLevel

    !> \brief Returns whether or not a message will be filtered
    procedure, private:: filter

    !> \brief Describes the object to its stream
    procedure, public:: describe
    
    ! >>> SETTERS

    !> \brief Sets the message filter used
    procedure, public:: setLevel

    !> \brief Sets the stream used
    procedure, public:: setStream

    !> \brief Returns a valid message level type
    procedure, private:: validateLevel

    !> \brief Resizes the filters used by the Logger
    procedure, private:: resizeFilters

    !> \brief Returns a flag stating if the filter exists
    procedure, private:: filterExists

    !> \brief Adds a vector of filters to the object
    procedure, public:: addFilters

    !> \brief Adds a single filter to the object
    procedure, public:: addFilter
    

    ! >>> GENERIC PROCEDURES

    !> \brief Returns the filter used for messages
    generic, public:: level => getLevel

    !> \brief Returns the stream unit
    generic, public:: stream => getStream

    !> \brief Returns the number of filters present
    generic, public:: numFilters => getNumFilters

    !> \brief Returns the the filters used
    generic, public:: filters => getFilters

end type Logger

