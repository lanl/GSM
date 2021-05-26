! =============================================================================
!>
!> \file
!> \brief  Contains the "buildMessage<>" procedures for derived Exceptions
!> \author CMJ
!>
! =============================================================================
  
  
  ! ===========================================================================
  !> \brief Build the standard exception's message
  subroutine buildMessage(this, errText)
    implicit none
    class(FException), intent(inout) :: this
    character(*),      intent(in   ) :: errText
    
    this%errorText = fmtGeneral // errText
    
  end subroutine buildMessage
  
  
  ! ===========================================================================
  !> \brief Build the unreachable exception's message
  subroutine buildMessageUnReachable(this, errText)
    implicit none
    class(UnReachable), intent(inout) :: this
    character(*),       intent(in   ) :: errText
    
    this%errorText = fmtUnReachable // errText
    
  end subroutine buildMessageUnReachable
  
  
  ! ===========================================================================
  ! Build the not implemented exception's message
  subroutine buildMessageUnImplemented(this, errText)
    implicit none
    class(NotImplemented), intent(inout) :: this
    character(*),          intent(in   ) :: errText
    
    this%errorText = fmtNotImplemented // errText
    
  end subroutine buildMessageUnImplemented
  
  
  ! ===========================================================================
  ! Build the bad allocation exception's message
  subroutine buildMessageBadAllocation(this, errText)
    implicit none
    class(BadAllocation), intent(inout) :: this
    character(*),         intent(in   ) :: errText
    
    this%errorText = fmtBadAllocation // errText
    
  end subroutine buildMessageBadAllocation
  
  
  ! ===========================================================================
  ! Build the under flow exception's message
  subroutine buildMessageUnderFlow(this, errText)
    implicit none
    class(UnderFlowError), intent(inout) :: this
    character(*),          intent(in   ) :: errText
    
    this%errorText = fmtUnderFlow // errText
    
  end subroutine buildMessageUnderFlow
  
  
  ! ===========================================================================
  ! Build the over flow exception's message
  subroutine buildMessageOverFlow(this, errText)
    implicit none
    class(OverFlowError), intent(inout) :: this
    character(*),         intent(in   ) :: errText
    
    this%errorText = fmtOverFlow // errText
    
  end subroutine buildMessageOverFlow
  
  
  ! ===========================================================================
  ! Build the logic error exception's message
  subroutine buildMessageLogicError(this, errText)
    implicit none
    class(LogicError), intent(inout) :: this
    character(*),      intent(in   ) :: errText
    
    this%errorText = fmtLogicException // errText
    
  end subroutine buildMessageLogicError
  
  
  ! ===========================================================================
  ! Build the InvalidArgument exception's message
  subroutine buildMessageInvalidArgument(this, errText)
    implicit none
    class(InvalidArgument), intent(inout) :: this
    character(*),           intent(in   ) :: errText
    
    this%errorText = fmtInvalidArg // errText
    
  end subroutine buildMessageInvalidArgument
  
  
  ! ===========================================================================
  ! Build the OutOfRange exception's message
  subroutine buildMessageOutOfRange(this, errText)
    implicit none
    class(OutOfRange), intent(inout) :: this
    character(*),      intent(in   ) :: errText
    
    this%errorText = fmtOutofRange // errText
    
  end subroutine buildMessageOutOfRange
  
  
  ! ===========================================================================
  ! Build the RunTime exception's message
  subroutine buildMessageRunTime(this, errText)
    implicit none
    class(RunTimeError), intent(inout) :: this
    character(*),        intent(in   ) :: errText
    
    this%errorText = fmtRunTime // errText
    
  end subroutine buildMessageRunTime
  
