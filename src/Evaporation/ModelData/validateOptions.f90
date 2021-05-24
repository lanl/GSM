
  subroutine validateOptions(evapData)

! ======================================================================
!
! Validates the options of the evaporation data class
!
!
! Written by CMJ, XCP-3, 8/2018 (Evap Class creation)
!
! ======================================================================

    use evaporationFissionData, only: levelDensityGCCIFlag, &
         & levelDensitySpecFlag

    implicit none
    class(EvaporationData), intent(inout) :: evapData

! ======================================================================

    if ( evapData%options%alev /= levelDensityGCCIFlag .and. &
         & evapData%options%alev < levelDensitySpecFlag ) then
       ! Use default value
       write(evapData%io%message, 1000) evapData%options%alev
       call evapData%io%print(2, 3, evapData%io%message)
       evapData%options%alev = defaultAlev
       write(evapData%io%message, 1010) evapData%options%alev
       call evapData%io%print(2, 3, evapData%io%message)
    end if


    return
! ======================================================================
1000 format("An invalid level density parameterization flag was ", &
          & "specified (", f6.2, ").")
1010 format("   Using the default parameterization flag of ", f6.2, ".")
! ======================================================================
  end subroutine validateOptions
