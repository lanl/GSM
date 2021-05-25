
  subroutine swapNuclei(gsmObj, proj, targ)

! ====================================================================
!
! Swaps a projectile and target nucleus' information
!
!
! Written by CMJ, XCP-3 (04/2019)
!
! ====================================================================

    implicit none
    class(GSM),            intent(inout) :: gsmObj
    class(GSMProjectile),  intent(inout) :: proj
    class(GSMTarget),      intent(inout) :: targ

    real(real64)      :: temp
    character(len= 4) :: pname
    character(len=10) :: name

! ====================================================================

    ! Write that the swap is occurring
    write(gsmObj%io%message, 1000) proj%numBaryons, proj%numProtons, &
            & nint(targ%numBaryons), nint(targ%numProtons)
    call gsmObj%io%print(3, 4, gsmObj%io%message)
    write(gsmObj%io%message, 1050)
    call gsmObj%io%print(3, 4, gsmObj%io%message)

    ! Swap values:
    ! (name)
    name = trim(adjustl(targ%particleName))
    targ%particleName = trim(adjustl(proj%particleName))
    proj%particleName = trim(adjustl(name))
    ! (A)
    temp = targ%numBaryons
    targ%numBaryons = dble( proj%numBaryons )
    proj%numBaryons = nint( temp, int32 )
    ! (Z)
    temp = targ%numProtons
    targ%numProtons = dble( proj%numProtons )
    proj%numProtons = nint( temp, int32 )
    ! (rest mass)
    temp = targ%restMass
    targ%restMass = proj%restMass
    proj%restMass = temp
    ! (A_f)
    temp = targ%afMultiplier
    targ%afMultiplier = proj%afMultiplier
    proj%afMultiplier = temp
    ! (C_Z)
    temp = targ%czMultiplier
    targ%czMultiplier = proj%czMultiplier
    proj%czMultiplier = temp

    ! Reverse the system:
    if ( proj%system == labSystem ) then
       proj%system = antiLabSystem
    else if ( proj%system == antiLabSystem ) then
       proj%system = labSystem
    else
       ! Shouldn't get here - revert to default system
       proj%system = defaultSystem
    end if

    return
! ====================================================================
1000 format("The projectile (A=", i3, ", Z=", i3, ") is larger than ", &
          & "the target (A=", i3, ", Z=", i3, ")")
1050 format("   The anti-lab. system will be used.")
! ====================================================================
  end subroutine swapNuclei
