
  function new_DecayChannels (modelDecayIN) result(obj)

! ======================================================================
!
! Constructs a new "DecayChannels" object
!
! ======================================================================
    implicit none
    logical,              intent(in   ), optional :: modelDecayIN
    type(decayChannels) :: obj

    logical :: modelDecay = .true.
    integer(int32) :: allocateStatus, i, j

! ======================================================================

    ! Set defaults
    if (present(modelDecayIN)) modelDecay = modelDecayIN

    ! Set and allocate
    obj%modelDecayOpt = modelDecay
    call obj%allocateData()

    ! If decay will be modelled, copy default values
    if (obj%modelDecay()) then
        ! Copy lookup mappings
        lookup: do i = 1, size(obj%lookDat, 1)
           obj%lookDat(i) = defaultChannels%look(i)
        end do lookup

        channelProbability: do i = 1, size(obj%cbrDat, 1)
           obj%cbrDat(i) = defaultChannels%cbr(i)
        end do channelProbability

        modesJ: do j = 1, size(obj%modeDat, 2)
           modesI: do i = 1, size(obj%modeDat, 1)
              obj%modeDat(i, j) = defaultChannels%mode(i, j)
           end do modesI
        end do modesJ

        ! TODO: Apply forced channels
    end if
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



