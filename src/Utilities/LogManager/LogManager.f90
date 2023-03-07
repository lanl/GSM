! =============================================================================
!
!> \file
!> \brief  Contains the implementation of the LogManager object
!> \author CMJ (XCP-3; LANL)
!
! =============================================================================
!
! DATA TYPE DESCRIPTION
!
!> \class LogManager
!> \brief Defines the base LogManager object
!
!> Defines the base LogManager object used by GSM and its sub-models and
!> model-data to manage a set of Logger objects
!
! =============================================================================
type, public:: LogManager
    private

    !> \brief The array of Logger objects
    type(Logger), private, dimension(:), allocatable:: d_loggers

    !> \brief The number of Logger objects contained in the array
    integer(int32), private:: d_numLoggers = 0_int32

    !> \brief The size of the Logger array
    integer(int32), private:: d_maxLoggers = 0_int32

    !> \brief A string to be printed
    character(:), public, allocatable:: message

 contains
    private

    ! >>> LOGGING PROCEDURES

    !> @{
    !> \brief Base message logging procedures
    procedure, private:: criticalBase
    procedure, private:: errorBase
    procedure, private:: warningBase
    procedure, private:: infoBase
    procedure, private:: debugBase
    procedure, private:: printBase
    !> @}

    !> @{
    !> \brief Interfaces to the base printing procedures
    procedure, private:: criticalExtended
    procedure, private:: errorExtended
    procedure, private:: warningExtended
    procedure, private:: infoExtended
    procedure, private:: debugExtended
    procedure, private:: printExtended
    !> @}

    !> @{
    !> \brief Client interface procedures to log specific messages
    generic, public:: critical => criticalBase, criticalExtended
    generic, public:: fatal    => criticalBase, criticalExtended
    generic, public:: error    => errorBase,    errorExtended
    generic, public:: warning  => warningBase,  warningExtended
    generic, public:: warn     => warningBase,  warningExtended
    generic, public:: info     => infoBase,     infoExtended
    generic, public:: log      => infoBase,     infoExtended
    generic, public:: debug    => debugBase,    debugExtended
    generic, public:: print    => printBase,    printExtended
    !> @}


    ! >>> INTROSPECTION

    !> \brief Returns the number of loggers present
    procedure, public:: getNumLoggers

    !> \brief Returns the size of the d_loggers array
    procedure, public:: getMaxLoggers

    !> \brief Returns the ith Logger of the manager
    procedure, public:: getLogger

    !> \brief Describes the loggers contained here
    procedure, public:: describe

    
    ! >>> SETTERS

    !> \brief Resizes the d_logger array
    procedure, public:: resize

    !> \brief Adds a logger to the array
    procedure, public:: addLogger
    generic,   public:: addLog => addLogger
    

    ! >>> GENERIC PROCEDURES

    !> \brief Returns the number of loggers present in the manager
    generic, public:: numLoggers => getNumLoggers

    !> \brief Returns the size of the d_loggers array
    generic, public:: maxLoggers => getMaxLoggers

    !> \brief Returns the i^th Logger object
    generic, public:: logger => getLogger

end type LogManager

