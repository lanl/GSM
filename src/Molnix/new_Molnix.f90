
  function new_Molnix(clientOptions, clientIO) result(molnixObj)

! =====================================================================
!
! Constructor for the Molnix class
!
! USE:
!    myMolnixObject = Molnix([newCevap], [newIOunit])
!
! REQUIRED ARGUMENTS: There are NO required arguments for this class
!
! OPTIONAL ARGUMENTS:
! (1) Use of 'newCevap' may be specified. The default value is 12.0.
!     'cevap' is used in the Molnix class as a multiplier. If the user
!     attempts to provide a value less than 0, the default will be used.
! (2) Use of 'iounit' may be specified. The default is STD_OUT (from the
!     intrinsic fortran module, iso_fortran_env). Users may direct where
!     output goes. Output is currently restricted to math warnings (namely
!     overflow errors) and comments.
!
!
! Written by CMJ, XCP-3, August 2018 (Molnix Class creation)
!
! =====================================================================

    implicit none
    type(molnixOptions),  intent(in   ), optional          :: clientOptions
    procedure(IOHANDLER), intent(in   ), optional, pointer :: clientIO
    type(Molnix)                                           :: molnixObj

! =====================================================================

    molnixObj%constructed = .TRUE.


    ! Set the I/O unit (if present)
    if ( present(clientIO) ) then
       if ( associated(clientIO) ) then
          molnixObj%io%print => clientIO
       else
          write(molnixObj%io%message, 1100)
          call molnixObj%io%print(2, 3, molnixObj%io%message)
       end if
    end if


    ! Set 'cevap' value (if present)
    if ( present(clientOptions) ) then

       molnixObj%options = clientOptions

       ! Verify valid options are used:
       if ( molnixOBj%options%cevap <=0 ) then
          molnixObj%options%cevap = defaultCevap
          write(molnixObj%io%message, 1000) clientOptions%cevap
          call molnixObj%io%print(2, 3, molnixObj%io%message)
          write(molnixObj%io%message, 1010) molnixObj%options%cevap
          call molnixObj%io%print(2, 3, molnixObj%io%message)
       end if
    end if


    return
! =====================================================================
1000 format("An invalid multiplier (", f8.3, ") was passed in to ", &
          & "the Molnix object.")
1010 format("   The default multiplier (", f8.3, ") will ", &
          & "be used instead.")
1100 format("The I/O handling procedure handed to the Molnix object ", &
          & "is not associated and will not be used.")
! =====================================================================
  end function new_Molnix
