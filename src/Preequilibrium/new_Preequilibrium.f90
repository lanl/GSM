
  function new_Preequilibrium(preeqData, clientRNG, clientFBU, clientOptions, &
       & clientPhotonEmission, clientIO) result(preeqObj)

! ======================================================================
!
! Constructor for the "Preequilibrium" class
!
! USE:
!    preequilibriumObject = Preequilibrium(fragmentBnk, preeqData, &
!         & newRNG, [npreqtyp], [r0], [ifam], [sigpre], [noprec], [preeqIO])
!
!
! REQUIRED ARGUMENTS: The fragment bank for progeny storage must be specified,
!    in addition to a Preequilibrium Data class and a RNG procedure pointer.
!
! OPTIONAL ARGUMENTS:
! (1) Npreqtyp may be specified to limit the consideration of the number
!     of particles that may be emitted from the residual. The default is
!     the maximum value (66) at present.
! (2) R0 may be specified to set a radius parameter that is used by
!     the class. Default value is 1.2d0.
! (3) IFam may be specified to set the coefficients a the level density
!     parameter calculation
! (4) sigpre may be specified to set how quickly the exponential for preeq.
!     emission decays (i.e. reduces or increase its width)
! (5) noprec may be specified to ignore all preequilibrium physics. This
!     option ought to never be used when simulating reality. A warning is
!     printed when this option is used.
! (6) PreeqIO may be specified so the client program has control of all
!     I/O that occurs during the Preequilibrium simulation. Currently, I/O
!     is restricted to comments, warnings, and errors. The I/O procedure
!     pointer that is allowed to be used has the arguments:
!
!        subroutine preeqIO ( [verbosity], [messageType], [message])
!
!     Where the [verbosity] indicates the "importance" of the message,
!     the [messageType] indicates the kind of message (comment, warning,
!     error, or fatal error), and the [message] is a string which contains
!     the actual message.
!
!
! Written by CMJ, XCP-3, 8/2018 (Class creation)
!
! ======================================================================

    use fermiBreakUpClass,       only: FermiBreakup, newFermiBreakup, &
         & fermiBreakUpOptions
    use preequilibriumData,      only: preeqDataDeclared
    use preequilibriumDataClass, only: PreequilibriumData

    implicit none
    type(PreequilibriumData),     intent(in   ), target   :: preeqData
    procedure(RANDOM),            intent(in   ), pointer  :: clientRNG
    class(FermiBreakup),          intent(in   ), optional :: clientFBU
    type(preequilibriumOptions),  intent(in   ), optional :: clientOptions
    procedure(PHOTOEMISSION),     intent(in   ), optional, pointer :: clientPhotonEmission
    procedure(IOHANDLER),         intent(in   ), optional, pointer :: clientIO
    type(Preequilibrium)                                  :: preeqObj

    type(fermiBreakUpOptions) :: fbuOptions

! ======================================================================

    ! Flag that object was successfully created:
    preeqObj%constructed = .TRUE.


    ! Determine if preequilibrium data has been established
    if ( .not. preeqDataDeclared ) then
       call preequilibriumInit()
    end if


    ! Point to the appropriate data:
    preeqObj%preeqData => preeqData
    if ( .not.preeqObj%preeqData%properlyConstructed() ) then
       write(preeqObj%io%message, 1300)
       call preeqObj%io%print(1, 2, preeqObj%io%message)
       preeqObj%constructed = .FALSE.
    end if
    ! (obtain local copies of the data)
    preeqObj%molEnergy => preeqObj%preeqData%getMolnix()
    preeqObj%fissBarr  => preeqObj%preeqData%getFissionBarrier()


    ! Point to the client's random number generator:
    if ( associated(clientRNG) ) then
       preeqObj%rng => clientRNG
    else
       write(preeqObj%io%message, 1100)
       call preeqObj%io%print(1, 2, preeqObj%io%message)
       preeqObj%constructed = .FALSE.
    end if


    ! Set where all printing goes (if client-defined):
    if ( present(clientIO) ) then
       if ( associated(clientIO) ) then
          preeqObj%io%print => clientIO
       else
          write(preeqObj%io%message, 1200)
          call preeqObj%io%print(2, 3, preeqObj%io%message)
       end if
    end if


    ! Set Fermi Break up object:
    if ( present(clientFBU) ) then
       preeqObj%fbuObj = clientFBU
       if ( .not. preeqObj%fbuObj%properlyConstructed() ) then
          preeqObj%fbuObj = newFermiBreakUp( preeqObj%rng, clientFBU%queryOptions(), &
               & preeqObj%io%print )
       end if
    else
       preeqObj%fbuObj = newFermiBreakUp( preeqObj%rng, fbuOptions, &
            & preeqObj%io%print )
    end if


    ! Set photon emission (if client-defined):
    if ( present(clientPhotonEmission) ) then
       if ( associated(clientPhotonEmission) ) then
          preeqObj%usePhotonEmission =  .TRUE.
          preeqObj%photonEmission    => clientPhotonEmission
       else
          write(preeqObj%io%message, 1000)
          call preeqObj%io%print(2, 3, preeqObj%io%message)
       end if
    end if


    ! Set options if desired by the client:
    if ( present(clientOptions) ) then

       ! Use client's options and validate them:
       preeqObj%options = clientOptions
       call preeqObj%validateOptions()

    end if


    return
! ======================================================================
1000 format("The photon emission procedure for the Preequilibrium object", &
          & "is not associated and will not be used.")
1100 format("The RNG procedure pointer given to the Preeq. object is ", &
          & "not associated and cannot be used.")
1200 format("The I/O procedure pointer given to the Preeq. object is ", &
          & "not associated and cannot be used.")
1300 format("The preeq. data object was not constructed at the ", &
          & "time of preeq. object construction.")
! ======================================================================
  end function new_Preequilibrium
