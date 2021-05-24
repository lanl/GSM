
  subroutine validateOptions(fbuObj)

! ====================================================================
!
! This procedure validates the options that clients have specified
! for the FBU object to use (i.e. it ensures their validity)
!
!
! Written by CMJ, XCP-3 (03/2019)
!
! ====================================================================

    implicit none
    class(FermiBreakUp), intent(inout) :: fbuObj

! ====================================================================

    ! No checks needed for akpScaling:
    if ( fbuObj%options%akpScalingFlag > 0 ) then
       write(fbuObj%io%message, 2000) fbuObj%options%akpScalingFlag
       call fbuObj%io%print(4, 4, fbuObj%io%message)
    end if


    ! Check recommended number of nucleons:
    if ( fbuObj%options%recNumNucleons > maxAllowedNucleons .or. &
         & fbuObj%options%recNumNucleons < minAllowedNucleons ) then
       write(fbuObj%io%message, 1000) fbuObj%options%recNumNucleons, &
            & minAllowedNucleons, maxAllowedNucleons
       call fbuObj%io%print(2, 3, fbuObj%io%message)
       fbuObj%options%recNumNucleons = defaultRecNumNucleons
       write(fbuObj%io%message, 1010) fbuObj%options%recNumNucleons
       call fbuObj%io%print(2, 3, fbuObj%io%message)
    end if



    return
! ====================================================================
1000 format("The requested nucleus size for Fermi Break-up (A=", i2, &
          & ") is outside the valid range [", i2, ", ", i2, "].")
1010 format("   Using default nucleus size (A=", i2, ").")
2000 format("Scaling of 'akp' (=", f7.3, ") will be utilized.")
! ====================================================================
  end subroutine validateOptions
