
  function new_PreequilibriumResults(fragmentBnk) result(results)

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

    implicit none
    type(preequilibriumFragment), intent(in   ), dimension(:), target :: fragmentBnk
    type(preequilibriumResults)   :: results

! ======================================================================

    ! Point to progeny bank
    results%progenyBnk => fragmentBnk


    ! Obtain size of progeny bank:
    results%maxProgeny = size(results%progenyBnk)


    ! Flag that object was successfully created:
    results%constructed = .TRUE.


    return
! ======================================================================
  end function new_PreequilibriumResults
