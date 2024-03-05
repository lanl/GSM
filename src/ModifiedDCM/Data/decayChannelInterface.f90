
  function new_DecayChannels (modelDecayIN, &
                              & forcedChannelIDs, &
                              & forcedChannelDecays) result(obj)

! ======================================================================
!
! Constructs a new "DecayChannels" object
!
! ======================================================================
    implicit none
    logical,        intent(in), optional :: modelDecayIN
    integer(int32), intent(in), optional, dimension(:) :: forcedChannelIDs
    integer(int32), intent(in), optional, dimension(:,:) :: forcedChannelDecays
    type(decayChannels) :: obj

    logical :: modelDecay = .true.

    ! Variables to aid in interim logic
    integer(int32), parameter :: maxLoop = 60
    integer(int32) :: &
        & numForcedChannels, &   ! Num. forced channels provided
        & numForcedModes, &      ! Num. forced decay modes provided
        & suggestedLoop         ! Suggested index for forced channel mapping
    ! Dummy variables only needed for other calls or looping
    integer(int32) :: i, j, ifl1, ifl2, ifl3, jspin, indx

! ======================================================================
    integer(int32) :: nforce, iforce, modes
    common/force/nforce,iforce(20),modes(5,20)
! ======================================================================

    ! Set defaults
    if (present(modelDecayIN)) modelDecay = modelDecayIN

    ! Set and allocate
    obj%modelDecayOpt = modelDecay
    call obj%allocateData()

    ! If decay will be modelled, copy default values
    useDecay: if (obj%modelDecay()) then
        ! Copy lookup mappings via interface
        lookup: do i = 1, size(obj%lookDat, 1)
           obj%lookDat(i) = defaultChannels%look(i)
        end do lookup

        ! Copy decay channel IDs via interface
        channelProbability: do i = 1, size(obj%cbrDat, 1)
           obj%cbrDat(i) = defaultChannels%cbr(i)
        end do channelProbability

        ! Copy channel progeny IDs via interface
        modesJ: do j = 1, size(obj%modeDat, 2)
           modesI: do i = 1, size(obj%modeDat, 1)
              obj%modeDat(i, j) = defaultChannels%mode(i, j)
           end do modesI
        end do modesJ
        obj%numEntriesDat = defaultChannels%numEntriesDat

        ! Overwrite data with specified forced channels
        forcedChannelsProvided: if (present(forcedChannelIDs) .AND. &
            & present(forcedChannelDecays)) then

           ! Check if provided arrays are of ok size / shape
           numForcedChannels = MIN(size(forcedChannelIDs), size(forcedChannelDecays, 2))
           numForcedModes = MIN(size(forcedChannelDecays, 1), 5)
           if (numForcedChannels <= 0 .OR. numForcedModes <= 0) then
              write(dataIO%message, "(A)") "Forced channels will not be applied, array sizes invalid."
              call dataIO%print(2, 3, dataIO%message)
              exit forcedChannelsProvided
           end if

           ! Set forced decay modes, if any specified and array size remains
           mappingsAvailable: if(obj%numEntriesDat < maxLoop) then
              forcedModes: do i = 1, numForcedChannels
                 obj%numEntriesDat = obj%numEntriesDat + 1_int32

                 ! Check if capacity reached; decrement and exit if so
                 forcedExceeded: if (obj%numEntriesDat > maxLoop) then
                    obj%numEntriesDat = obj%numEntriesDat - 1_int32
                    exit forcedModes
                 else
                    ! Determine the particle indx for the associated particle ID who's
                    ! channel has a forced decay specified
                    call flavor(forcedChannelIDs(i), ifl1, ifl2, ifl3, jspin, indx)

                    ! TODO: Instead of growing the mapping tables, we ought to obtain
                    ! the existing mapping, if one exists, and then override the mode values.
                    ! with cbr = 1, any sampling of that particle's decay will force
                    ! usage of that decay channel.
                    !
                    ! This would look like the following unreachable code
                    ! ("not_implemented" contract added for safety). This was not
                    ! implemented since it could not be easily tested and verified
                    ! against real data at the time of consideration
                    if (.false.) then
                       call not_implemented("Map restriction", &
                          & __FILE__, &
                          & __LINE__)
                       suggestedLoop = obj%look(indx)
                       if (suggestedLoop == 0) then
                          ! ID is unmapped; create a new mapping
                          ! (otherwise, use existing)
                          suggestedLoop = obj%numEntriesDat
                          obj%lookDat(indx) = suggestedLoop
                       end if
                    else
                       ! Add particle ID lookup table and associated values
                       suggestedLoop = obj%numEntriesDat
                       obj%lookDat(indx) = suggestedLoop
                    end if

                    ! Now force the channel and modes
                    ! NOTE: "cbr" is setup as a channel probability associated with
                    !       a random sampling; setting = 1 will force that channel.
                    !       If we overwrite a mapping, subsequent data rows for
                    !       that channel will be abandoned, if any exist.
                    obj%cbrDat(suggestedLoop) = one
                    obj%modeDat(1 : numForcedModes, suggestedLoop) = forcedChannelDecays(1 : numForcedModes, i)
                    if (numForcedModes < 5) then
                       obj%modeDat(numForcedModes : 5, suggestedLoop) = 0_int32
                    end if
                 end if forcedExceeded
              end do forcedModes
           end if mappingsAvailable
        end if forcedChannelsProvided
    end if useDecay
    return
  end function new_DecayChannels

  subroutine allocateData(obj)

! ======================================================================
!
! Allocates and initializes the internal arrays of a decayChannel object.
!
! ======================================================================

    implicit none
    class(decayChannels), intent(inout) :: obj

    integer(int32) :: state = 0_int32

! ======================================================================

    ! Allocate lookup table
    if(allocated(obj%lookDat)) deallocate(obj%lookDat)
    allocate(obj%lookDat(400), source=0_int32, stat=state)
    Insist(state == 0, "Lookup table allocation failed")

    ! Allocate probability table
    if(allocated(obj%cbrDat)) deallocate(obj%cbrDat)
    allocate(obj%cbrDat(600), source=1.0_real64, stat=state)
    Insist(state == 0, "Channel probability allocation failed")

    ! Allocate decay progeny ID table
    if(allocated(obj%modeDat)) deallocate(obj%modeDat)
    allocate(obj%modeDat(5, 600), source=0_int32, stat=state)
    Insist(state == 0, "Decay channel progeny list allocation failed")

    return
  end subroutine allocateData


  function modelDecay (obj) result(dataValue)

! ======================================================================
!
! Interface for the "modelDecay" class option; this will allow consumers
! to evaluate if they are modelling decay or not
!
! ======================================================================
    implicit none
    class(decayChannels), intent(in   ) :: obj
    logical :: dataValue

    dataValue = obj%modelDecayOpt
    return
  end function modelDecay


  function look (obj, i) result(dataValue)

! ======================================================================
!
! Interface for the "lookDat" variable; this is the lookup mapping table
! that maps a particle ID to the data utilized
!
! ======================================================================

    implicit none
    class(decayChannels), intent(in   ) :: obj
    integer(int32),       intent(in   ) :: i
    integer(int32) :: dataValue

    Require(0 < i .AND. i <= size(obj%lookDat, 1))
    dataValue = obj%lookDat(i)
    return
  end function look


  function cbr (obj, i) result(dataValue)

! ======================================================================
!
! Interface for the "cbrDat" variable; this is effectively the cumulative
! decay channel probability when randomly sampled.
!
! This is used to randomly sample the decay channel.
!
! ======================================================================

    implicit none
    class(decayChannels), intent(in   ) :: obj
    integer(int32),       intent(in   ) :: i
    real(real64) :: dataValue

    Require(0 < i .AND. i <= size(obj%cbrDat, 1))
    dataValue = obj%cbrDat(i)
    return
  end function cbr


  function mode (obj, i, j) result(dataValue)

! ======================================================================
!
! Interface for the "modeDat" variable; this indicates the secondary
! particles emitted from the decay of the parent. mode values are the ID
! of the secondary, which can additionally be mapped similarly again
! for subsequent decays.
!
! ======================================================================

    implicit none
    class(decayChannels), intent(in   ) :: obj
    integer(int32),       intent(in   ) :: i
    integer(int32),       intent(in   ) :: j
    real(real64) :: dataValue

    Require(0 < i .AND. i <=   size(obj%modeDat, 1))
    Require(0 < j .AND. j <=   size(obj%modeDat, 2))
    dataValue = obj%modeDat(i, j)
    return
  end function mode



