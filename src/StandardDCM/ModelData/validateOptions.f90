
  subroutine validateOptions(dataObj)

! ==============================================================================
!
! Validates all options used by the SDCMData object
!
!
! Written by CMJ, XCP-3, 12/2018
!
! ==============================================================================

    use, intrinsic :: iso_fortran_env, only: int32
    use standardDCMDataParams, only: zro, one

    implicit none
    class(standardDCMData), intent(inout) :: dataObj

    integer(int32) :: i

! ==============================================================================

    ! Verify only physical values are present
    ! (number of zones)
    if ( dataObj%options%numZones <= minZonesAllowed .or. &
         & dataObj%options%numZones >= maxZonesAllowed ) then
       write(dataObj%io%message, 1100) dataObj%options%numZones
       call dataObj%io%print(0, 3, dataObj%io%message)
       write(dataObj%io%message, 1110) minZonesAllowed, maxZonesAllowed
       call dataObj%io%print(0, 3, dataObj%io%message)
       write(dataObj%io%message, 1120) numZonesDefault
       call dataObj%io%print(0, 3, dataObj%io%message)
       dataObj%options%numZones = numZonesDefault
    end if

    ! (exponential term)
    if ( dataObj%options%r0 <= zro ) then
       write(dataObj%io%message, 1200) dataObj%options%r0
       call dataObj%io%print(0, 3, dataObj%io%message)
       write(dataObj%io%message, 1999) r0Default
       call dataObj%io%print(0, 3, dataObj%io%message)
       dataObj%options%r0 = r0Default
    end if
    if ( dataObj%options%maxRad <= dataObj%options%r0 ) then
       write(dataObj%io%message, 1300) dataObj%options%maxRad
       call dataObj%io%print(0, 3, dataObj%io%message)
       write(dataObj%io%message, 1999) maxRadDefault
       call dataObj%io%print(0, 3, dataObj%io%message)
       dataObj%options%maxRad = maxRadDefault
    end if
    if ( abs(dataObj%options%expDenom) <= 1.0d-20 ) then
       write(dataObj%io%message, 1400)
       call dataObj%io%print(0, 3, dataObj%io%message)
       write(dataObj%io%message, 1410) dataObj%options%expDenom
       call dataObj%io%print(0, 3, dataObj%io%message)
       write(dataObj%io%message, 1999) expDenomDefault
       call dataObj%io%print(0, 3, dataObj%io%message)
       dataObj%options%expDenom = expDenomDefault
    end if

    ! for the 'aveDen' term
    do i = 1, maxZonesAllowed
       if ( dataObj%options%aveDen(i) < zro .or. dataObj%options%aveDen(i) > one ) then
          ! Invalid value; warn client and change to default
          write(dataObj%io%message, 1500) dataObj%options%aveDen(i), i
          call dataObj%io%print(0, 3, dataObj%io%message)
          write(dataObj%io%message, 1510) aveDenDefault(i)
          call dataObj%io%print(0, 3, dataObj%io%message)
          dataObj%options%aveDen(i) = aveDenDefault(i)
       end if
    end do

    return
! ==============================================================================
1100 format("An invalid number of nuclear zones was requested (", i3, &
          & ").")
1110 format("   No less than ", i3, " and no more than ", i3, &
          & " zone(s) may be utilized.")
1120 format("   Defaulting to ", i3, " zone(s).")
1200 format("An invalid radius multipler was requested (", es11.4, &
          & " fm).")
1300 format("An unphysical maximum nuclear radius was selected (", es11.4, &
          & " fm).")
1400 format("The selected exponential denominator for the Standard DCM Data")
1410 format("   object cannot be used (", es11.4, " fm).")
1500 format("The chosen average nucleon density (", es11.4, ") for zone ", &
          & i3, " is not valid.")
1510 format("   Defaulting to ", es11.4, " for this zone.")
1999 format("   Defaulting to ", es11.4, " fm.")
! ==============================================================================
  end subroutine validateOptions
