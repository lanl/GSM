! =============================================================================
!>
!> \file
!> \brief  Contains the Exceptions types
!> \author CMJ
!>
! =============================================================================

  !> \brief The base Fortran Exception object
  !>
  !> The base Fortran Exception object that other Fortran Exception objects
  !> are derived off of.
  type, public :: FException
     private

     !> \brief The text produced by the exceiption (e.g. what it is and location)
     character(:), private, allocatable :: errorText
     
     !> \brief The location the exception was thrown in
     character(:), private, allocatable :: errorLocation
     
   contains
     private
     
     !> \brief Returns the location the exception was generated in
     procedure, public, non_overridable :: location
     
     !> \brief Returns the text of the exception
     procedure, public, non_overridable :: what
     
     !> \brief Builds the exception message
     procedure, private :: buildMessage
     
     !> \brief Creates the exception object
     procedure, public, non_overridable :: throw
     generic, public :: new => throw
     
  end type FException


  !> \brief Indicator of unreachable code being accessed
  type, public, extends(FException) :: UnReachable
     private
   contains
     private
     procedure :: buildMessage => buildMessageUnReachable
  end type UnReachable


  !> \brief Indicates a procedure has not been implemented
  type, public, extends(FException) :: NotImplemented
     private
   contains
     private
     procedure :: buildMessage => buildMessageUnImplemented
  end type NotImplemented


  !> \brief Indicates a bad allocation ocurred
  type, public, extends(FException) :: BadAllocation
     private
   contains
     private
     procedure :: buildMessage => buildMessageBadAllocation
  end type BadAllocation


  !> \brief Indicates an underflow error occurred
  type, public, extends(FException) :: UnderFlowError
     private
   contains
     private
     procedure :: buildMessage => buildMessageUnderFlow
  end type UnderFlowError


  !> \brief Indicates an overflow error occurred
  type, public, extends(FException) :: OverFlowError
     private
   contains
     private
     procedure :: buildMessage => buildMessageOverFlow
  end type OverFlowError


  !> \brief Indicates a general logic exception occurred
  type, public, extends(FException) :: LogicError
     private
   contains
     private
     procedure :: buildMessage => buildMessageLogicError
  end type LogicError
  

  !> \brief Indicates an invalid argument was received
  type, public, extends(FException) :: InvalidArgument
     private
   contains
     private
     procedure :: buildMessage => buildMessageInvalidArgument
  end type InvalidArgument


  !> \brief Indicates an Out of Range error occurred
  type, public, extends(FException) :: OutOfRange
     private
   contains
     private
     procedure :: buildMessage => buildMessageOutOfRange
  end type OutOfRange


  !> \brief Indicates a general run time error occurred
  type, public, extends(FException) :: RunTimeError
     private
   contains
     private
     procedure :: buildMessage => buildMessageRunTime
  end type RunTimeError

       
