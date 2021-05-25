! =============================================================================
!
!> \file
!> \brief  Contains the default values for the objects in the Loggers module
!> \author CMJ (XCP-3; LANL)
!
! =============================================================================

    !> \brief Defines the number of message types contained
    integer(int32), private, parameter:: numMessageLevels = 7_int32

    !> \brief Defines numerics for the logging filters/types
    enum, bind(c)
        enumerator:: QUIETLog     =  0_c_int
        enumerator:: CRITICALLog  = 10_c_int
        enumerator:: ERRORLog     = 20_c_int
        enumerator:: WARNINGLog   = 30_c_int
        enumerator:: INFOLog      = 40_c_int
        enumerator:: DEBUGLog     = 50_c_int
        enumerator:: UNSETLog     = 60_c_int
    end enum
    public:: QUIETLog, CRITICALLog, ERRORLog, WARNINGLog, INFOLog, &
        & DEBUGLog, UNSETLog

    !> @{
    !> Labels prepended to a message for different message types
    character(*), private, parameter:: preCritical = "(F)   Fatal: "
    character(*), private, parameter:: preError    = "(E)   Error: "
    character(*), private, parameter:: preWarning  = "(W) Warning: "
    character(*), private, parameter:: preInfo     = "       Info: "
    character(*), private, parameter:: preDebug    = "      Debug: "
    character(*), private, parameter:: preUnset    = ""
    !> }

    !> @{
    !> Formats used to print messages
    character(*), private, parameter:: stringFormat = "(A)"
    !> @}

    !> \brief Defines the default message level printed
    integer(c_int), private, parameter:: defaultLevel = UNSETLog

    !> \brief Defines the default stream messages are sent to
    integer(int32), private, parameter:: defaultStream = output_unit

