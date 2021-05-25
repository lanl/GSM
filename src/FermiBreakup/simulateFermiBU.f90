
  subroutine simulateFermiBU(fbuObj, residualNucleus, results)

! ====================================================================
!
! This procedure interfaces to the main FBU simulation routine (fermid)
!
!
! Written by CMJ, XCP-3 (03/2019)
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32
    use fermiBreakupParams, only: nucleon_mass

    implicit none
    class(FermiBreakUp),       intent(inout) :: fbuObj
    type(fermiBreakUpNucleus), intent(in   ) :: residualNucleus
    type(fermiBreakUpResults), intent(inout) :: results

    integer(int32) :: numBaryons = 0_int32

! ====================================================================

    ! Validate FBU object construction:
    if ( .not. fbuObj%constructed ) then
       write(fbuObj%io%message, 1000)
       call fbuObj%io%print(1, 2, fbuObj%io%message)
       write(fbuObj%io%message, 2000)
       call fbuObj%io%print(1, 2, fbuObj%io%message)
       results%simState = results%simState + fbuObjectNotConstructed
    end if

    ! For results object:
    if ( .not. results%constructed ) then
       write(fbuObj%io%message, 1100)
       call fbuObj%io%print(1, 2, fbuObj%io%message)
       write(fbuObj%io%message, 2000)
       call fbuObj%io%print(1, 2, fbuObj%io%message)
       results%simState = results%simState + fbuResultNotConstructed
    end if

    ! Check if array sizes *might* be exceeded based on nucleus size:
    numBaryons = nint( residualNucleus%numBaryons )
    if ( .not.fbuObj%recommendFermiBreakUp( numBaryons ) ) then
       write(fbuObj%io%message, 1200) fbuObj%options%recNumNucleons
       call fbuObj%io%print(1, 2, fbuObj%io%message)
       ! Check for allowance of client FBU usage:
       if ( minAllowedNucleons <= numBaryons .and. numBaryons <= maxAllowedNucleons ) then
          fbuObj%options%recNumNucleons = numBaryons
          write(fbuObj%io%message, 1210) fbuObj%options%recNumNucleons
          call fbuObj%io%print(1, 3, fbuObj%io%message)
       else
          ! Warn user:
          write(fbuObj%io%message, 1220) numBaryons
          call fbuObj%io%print(1, 2, fbuObj%io%message)
          ! Store nucleus as a ''progeny'':
          results%numProgeny = 1
          results%progenyBnk(1)%numBaryons = residualNucleus%numBaryons
          results%progenyBnk(1)%numProtons = residualNucleus%numProtons
          results%progenyBnk(1)%kinEnergy  = residualNucleus%kinEnergy
          results%progenyBnk(1)%restMass   = nucleon_mass * residualNucleus%numBaryons
          results%progenyBnk(1)%linearXMom = residualNucleus%linearXMom
          results%progenyBnk(1)%linearYMom = residualNucleus%linearYMom
          results%progenyBnk(1)%linearZMom = residualNucleus%linearZMom
          ! Flag the error:
          results%simState = results%simState + invalidNucleus
       end if
    end if

    ! If any of the above fatal warnings/errors are found exit:
    if ( results%simState /= noSimulationWarnings ) return


    ! Simulate FBU physics now:
    call fbuObj%fermid(residualNucleus, results)


    return
! ====================================================================
1000 format("The Fermi Break-Up object was not properly constructed.")
1100 format("The Fermi Break-Up results object was not properly constructed.")
1200 format("The Fermi Break-Up object is being used outside its ", &
          & "constructed specifications for the nucleus (A=", i3, ").")
1210 format("   The object will be expanded to encompass the ", &
          & "requested nucleus (A=", i3, ").")
1220 format("   The requested nucleus (A=", i3, ") CANNOT be ", &
          & "disintegrated via Fermi Break-up physics.")
2000 format("  Unable to simulate Fermi Break-Up physics.")
! ====================================================================
  end subroutine simulateFermiBU



