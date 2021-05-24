
  subroutine validateOptions ( evapObj )

! ===================================================================================
!
! This subroutine validates the options that the evaporation object is to use
! during its simulations.
!
!
! Written by CMJ, XCP-3, 01/2019
!
! ===================================================================================

    implicit none
    class(Evaporation), intent(inout) :: evapObj

! ===================================================================================


    ! Check that number of evaporated particles is within valid range
    ! RANGE: 1 <= numEvapType <= 66
    if ( evapObj%options%numEvapType < minNumEvapType .or. &
         & evapObj%options%numEvapType > maxNumEvapType ) then
       ! Value is outside the valid range; revert to default and warn user
       write(evapObj%io%message,1010) evapObj%options%numEvapType
       call evapObj%io%print( 2, 3, evapObj%io%message)

       evapObj%options%numEvapType = defaultNumEvapType

       write(evapObj%io%message,1000) evapObj%options%numEvapType
       call evapObj%io%print( 2, 3, evapObj%io%message)
    end if



    ! Set inverse reaction cross section parameter option
    ! Valid values:
    ! =  0         :: Precis Parameter Set
    ! 1 < # < 10   :: Simple set (r0 =   #)
    ! = 10         :: Simple set (r0 = 1.5)
    ! = -1         :: Furihata-Botniva Set
    ! <  0 (/= -1) :: Furihata-Atchison Set
    if ( evapObj%options%inverseParameter > maxInverseParameter ) then
       ! Outside valid range - warn client and correct
       write(evapObj%io%message, 1100) evapObj%options%inverseParameter
       call evapObj%io%print( 2, 3, evapObj%io%message)

       evapObj%options%inverseParameter = defaultInverseParameter

       write(evapObj%io%message, 1110) evapObj%options%inverseParameter
       call evapObj%io%print( 2, 3, evapObj%io%message)
    end if




    ! Validate fission parameterization used:
    ! Valid values:
    ! <= 0      :: Furihata's Updated Parameterization
    ! >  0      :: Original Atchison Set
    ! NOTE: NO check is needed.




    ! Validate the use of reduced computation time:
    ! Valid values:
    ! <= 0   :: Use standard simulation
    ! >  0   :: Use reduced computation time option
    ! NOTE: NO check is needed.


    return

! ===================================================================================
1000 format(i4, " particles cannot be considered for evaporation.")
1010 format("   ", i3, " will instead be considered (default).")
1100 format("The inverse cross section parameterization set flag is invalid (", &
          & f5.2, ").")
1110 format("   Utilizing the default parameterization (", f5.2, ").")
! ===================================================================================
  end subroutine validateOptions
