! =============================================================================
!
!> \file
!> \brief  Contains the tstLoggerModule, testing the LogManger object
!> \author CMJ (XCP-3; LANL)
!
!> Contains the tests for the LogManager object
!
! =============================================================================
module tstLogManager

    ! Import the Logger object and its constructor
    use Loggers, only: Logger, newLogger
    use LogManagers, only: LogManager, newLogManager

    implicit none
    private

    type(LogManager), public:: io

    public:: tstSizing
    public:: tstGetters
    public:: tstMessaging

 contains

     subroutine tstSizing()
         implicit none
         type(Logger) :: logType
         type(Logger), dimension(15) :: lotsLoggers

         ! Initialize with one logger
         io = newLogManager(lotsLoggers)

         ! Add a few loggers...
         call io%addLogger(logType)
         call io%addLogger(logType)

         ! Remove all loggers
         call io%resize(-1)
         
         ! Add some loggers again
         call io%addLogger(logType)
         call io%resize(10)
         call io%addLogger(logType)

         ! Empty manager now:
         call io%resize(0)
         
         return
     end subroutine tstSizing


     subroutine tstGetters()
         implicit none
         type(Logger):: newLogger

         ! Number of loggers
         write(*,"(A, i3)") "Number of loggers: ", io%getNumLoggers()

         ! Get a logger
         call io%resize(-1)
         newLogger = io%logger(3)

         return
     end subroutine tstGetters


     subroutine tstMessaging()
         implicit none
         type(Logger):: first
         type(Logger):: second
         type(Logger):: other

         ! Construct Loggers
         first = newLogger()
         second = newLogger()
         other = newLogger()

         ! Construct LogManager
         io = newLogManager()
         call io%addLog(first)
         call io%addLog(second)

         ! Now test
         call io%info("Here's a message!")

         io%message = &
             & "Note that it's form is bad usage as race conditions " // &
             & "(in parallel) may occur!"
         call io%debug()


         call io%print(&
             & "Use a message within the call if at all possible, please.")
         
         write(io%message, "(A)") "Testing the likely usage of clients..."
         call io%warning()
         
         return
     end subroutine tstMessaging


end module tstLogManager

