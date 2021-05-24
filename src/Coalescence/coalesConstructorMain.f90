
  function coalesConstructorMain( clientData, clientOptions, &
       & clientIO) result(coalObj)

! ======================================================================
!
! Constructor for Coalescence object.
!
! USE:
!    coalsObj = newCoalescence(data, [options], [newIO])
!
!
! REQUIRED ARGUMENTS: Constructor must be fed the arrays
! for both the coales_pmemo and coales_imemo data structures and the
! incident particles energy from the reaction.
!
! OPTIONAL ARGUMENTS:
! (1) Use of the expanded coalescence model (expCoal) can be specified,
!     where the default is to use the expanded coalescence model. The
!     (=0) Only coalesce up to d, t, He3, He4
!     (=2) Coalesce (=0) fragments, then He6, Li6, Li7, and Be7
! (2) The newIO allows client programs to specify and control the use of verbose I/O.
!     By default, the class contains its own I/O handler, however client programs
!     can control it via this procedure pointer.
! (3) The number of fragments to be considered for coalescence may be
!     included, limiting the number of produced fragments considered for
!     the coalescence process.
! (4) A list (dimension 4) of momentum spheres may be fed into the
!     constructor to define different momentum sphere values within the
!     class. The class will use the default value for that sphere
!     if the passed in value of the array is <0. Setting sphere to 0
!     will result in NO coalescence resulting in that particle (d/t,He3/He4/LF)
!
!
! Written by CMJ, XCP-3, July 2018
!
! ======================================================================

    implicit none
    type(CoalescenceData),    intent(in   ), optional :: clientData
    type(coalescenceOptions), intent(in   ), optional :: clientOptions
    procedure(IOHANDLER),     intent(in   ), optional, pointer :: clientIO
    type(Coalescence)                                 :: coalObj

! ======================================================================

    ! Determine where I/O goes
    if ( present(clientIO) ) then
       if ( associated(clientIO) ) then
          coalObj%io%print => clientIO
       else
          write(coalObj%io%message, 2000)
          call coalObj%io%print(2, 3, coalObj%io%message)
       end if
    end if


    ! Use client data object if defined:
    if ( present(clientData) ) then
       coalObj%data = clientData
       if ( .not. coalObj%data%properlyConstructed() ) then
          ! Client didn't construct data object - warn and continue
          write(coalObj%io%message, 1000)
          call coalObj%io%print(2, 3, coalObj%io%message)
          write(coalObj%io%message, 1010)
          call coalObj%io%print(2, 3, coalObj%io%message)
       end if
    end if


    ! Set use of expanded coalescence model
    if ( present(clientOptions) ) then
       ! (<0) signals NO use, (>=0) signals use
       coalObj%options = clientOptions
       ! Warn user (higher nucleon yields than is typically physical)
       if ( coalObj%options%expandedCoalescence <= 0 ) then
          write(coalObj%io%message,1200)
          call coalObj%io%print(3, 4, coalObj%io%message)
          write(coalObj%io%message,1210)
          call coalObj%io%print(3, 4, coalObj%io%message)
       end if
    end if


    ! Coalescence object was constructed successfully
    coalObj%constructed = .TRUE.

    return
! ======================================================================
1000 format("The coalescence data object was not constructed.")
1010 format("   Coalescence radii will be approximated as the defaults.")
1200 format("The expanded coalescence model will not be used.")
1210 format("   Larger nucleon yields will be expected.")
2000 format("The message handling procedure passed to the Coalescence ", &
          & "object is not associated and will not be used.")
! ======================================================================
  end function coalesConstructorMain
