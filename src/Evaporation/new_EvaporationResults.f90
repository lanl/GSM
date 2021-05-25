
  function new_EvaporationResults (fragmentBnk, fissionBnk) result(evapResults)

! ===================================================================================
!
! Constructor for the 'Evaporation' class.
!
! USE:
!    evaporationObject = Evaporation(fragmentBnk, evaporationData, [options], [clientIO] )
!
!
! REQUIRED ARGUMENTS: Constructor must be given the fragment bank array for the particle
! bank data structure (evaporationFragment) to store the fragment data in, AND data
! that used for the simulation, but that is based on a calculation option for the
! level density.
!
! OPTIONAL ARGUMENTS:
! (1) Clients may include the 'evaporationOptions' type if desired. Doing so allows
!     client programs to control how the evaporation model simulates and what
!     parameterizations to use.
! (2) The 'clientIO' variable allows clients to control how the evaporaiton model
!     outputs messages. The evaporation model makes calls to this procedure whenever
!     a message of any importance is created, giving the procedure a flag stating
!     how important the message is and a flag for the type of message (comment, error, etc.)
!
!
! Written by CMJ, XCP-3, 8/2018 (Class creation)
!
! ===================================================================================

    use evaporationDataClass,   only: EvaporationData

    implicit none
    type(evaporationFragment), intent(in   ), dimension(:), target :: fragmentBnk
    type(fissionFragment),     intent(in   ), dimension(:), target :: fissionBnk
    type(evaporationResults)                                       :: evapResults

! ===================================================================================

    ! Set progeny bank information:
    evapResults%progenyBnk        => fragmentBnk
    evapResults%progenyBnkSize    =  size ( evapResults%progenyBnk )
    evapResults%progenyBnkSizeM1  =  evapResults%ProgenyBnkSize - 1

    ! Set fission bank information:
    evapResults%fissionBnk     => fissionBnk
    evapResults%fissionBnkSize =  size ( evapResults%fissionBnk )



    evapResults%constructed = .TRUE.

    return
! ===================================================================================
  end function new_EvaporationResults
