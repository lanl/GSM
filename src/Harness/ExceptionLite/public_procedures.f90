!==============================================================================
!>
!> Contains procedure implementation for the ExceptionLite module.
!>
!==============================================================================

!==============================================================================
!
!> \fn getNewID
!> \brief Obtains a new ID for an exception
!>
!> Obtain a new ID for an exception object.
!
!==============================================================================
   function getNewID() result(id)
       ! Obtain the ID of the next exception object
       integer(int32) :: id = next_id

       ! Obtain a new ID to be used for the next exception
       ! (in this case, let's just be simple and increment)
       next_id = next_id + 1
   end function

!==============================================================================
!
!> \fn throwUnhandled
!> \brief Halt program execution (for when consumer doesn't handle an exception)
!>
!> Throws an exception if it is unhandled. This will update the global exception
!> stack as well to mark an exception as handled.
!>
!> ARGUMENTS:
!> \param[in] ex       An ExceptionLite class
!> \param[in] message  (optional) An additional string to print prior to exiting
!
!==============================================================================
   subroutine throwUnhandled(ex, message)
      class(ExceptionLite), intent(in) :: ex
      character(:), allocatable, optional, intent(in) :: message

      ! Write any additional message(s) provided by the client
      ! For example, if client wants to flush any stored I/O
      ! for delayed printing
      if (present(message)) then
          write(error_unit, "(A)") message
      end if

      ! Write the exception message
      write(error_unit, "A (A) exception occurred: (A)") ex%kind(), ex%what()
      
      ! Write location information if compiled to do so
      if (ExLocationDetails) then
          write(error_unit, " ^^^ Error at '(A)', line (i5):") ex%file, ex%line
      end if

      error stop "Execution terminating."
   end subroutine

!==============================================================================
!
!> \fn throwExceptionInternal
!> \brief Create an exceptions object, internal to the module
!>
!> Create an exception object that is internal to the ExceptionLite module.
!> This method should be called for all exception types that are thrown, again
!> being called by module method implementations (e.g., throwing a XXX
!> exception)
!>
!> ARGUMENTS:
!> \param[in] message  An additional string to print prior to exiting
!> \param[in] file     The file the exception occurred in
!> \param[in] line     The line the exception occurred on within the file
!> \param[in] exType   The exception type, based on the module enumeration.
!
!==============================================================================
 subroutine throwExceptionInternal(message, file, line, exType)
     character(:), allocatable, intent(in) :: message
     character(:), allocatable, intent(in) :: file
     integer(int16), intent(in) :: line
     integer(int16), intent(in) :: exType

     type(ExceptionLite) :: ex

     ! Validate parameters
     ! - Ensure the line is valid
     if (line <= 0) then
         call throwLogicException("Cannot use a line less than 0.", &
                  & __FILE__, &
                  & __LINE__)
     end if
     ! - Ensure the exception type is valid
     select case (ex%kind)
        case (GENERAL_EXCEPTION)
        case (LOGIC_EXCEPTION)
        case (RUNTIME_EXCEPTION)
        case default
            call throwLogicException("Invalid exception type, please review.", &
                & __FILE__, &
                & __LINE)
     end select

     id = getNewID()

     ! Parameters have been validated and appear correct; create an exception
     ! object using the defined constructor and add to the internal array.
     ex = ExceptionLite(id, message, file, line, exType)
     call addException(ex)
 end subroutine

!==============================================================================
!
!> \fn throwRunTimeException
!> \brief Creates a runtime exception
!>
!> Create a runtime exception using the specified message. This exception type
!> should be used for halting execution in the event of runtime issues, such
!> as memory leaks, mathemtical underflows and overflows, generic system events,
!> etc.
!>
!> ARGUMENTS:
!> \param[in] message  An additional string to print prior to exiting
!> \param[in] file     The file the exception occurred in
!> \param[in] line     The line the exception occurred on within the file
!
!==============================================================================
 subroutine throwRunTimeException(message, file, line)
     character(:), allocatable, intent(in) :: message
     character(:), allocatable, intent(in) :: file
     integer(int16), intent(in) :: line

     call throwExceptionInternal(message, file, line, RUNTIME_EXCEPTION)
 end subroutine

!==============================================================================
!
!> \fn throwLogicException
!> \brief Creates a logic exception
!>
!> Create a logic exception using the specified message. This exception type
!> should be used for halting execution in the event of failed program logic,
!> such as in parameter validation or physical constraints.
!>
!> ARGUMENTS:
!> \param[in] message  An additional string to print prior to exiting
!> \param[in] file     The file the exception occurred in
!> \param[in] line     The line the exception occurred on within the file
!
!==============================================================================
 subroutine throwLogicException(message, file, line)
     character(:), allocatable, intent(in) :: message
     character(:), allocatable, intent(in) :: file
     integer(int16), intent(in) :: line

     call throwExceptionInternal(message, file, line, LOGIC_EXCEPTION)
 end subroutine

!==============================================================================
!
!> \fn throwException
!> \brief Creates a basic exception object
!>
!> Create an exception using the specified message. This exception type
!> should be used for halting execution for program failures that do not fall
!> into another pre-defined category.
!>
!> ARGUMENTS:
!> \param[in] message  An additional string to print prior to exiting
!> \param[in] file     The file the exception occurred in
!> \param[in] line     The line the exception occurred on within the file
!
!==============================================================================
 subroutine throwException(message, file, line)
     character(:), allocatable, intent(in) :: message
     character(:), allocatable, intent(in) :: file
     integer(int16), intent(in) :: line

     call throwExceptionInternal(message, file, line, GENERAL_EXCEPTION)
 end subroutine

   !> Retrieve an exception
   public :: retrieveException

   !> Clear a single exception given a retrieved_exception
   public :: handleException

   !> Returns the ID of the next exception type and predicts the next exception
   !> ID
   private :: retrieveID

