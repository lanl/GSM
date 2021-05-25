
  subroutine validateOptions( preeqObj )

! ======================================================================
!
! Validates the options that the preequilibrium class is to use:
!
! ======================================================================

    implicit none
    class(Preequilibrium), intent(inout) :: preeqObj

! ======================================================================


    ! Verify radius parameter is ok:
    if ( preeqObj%options%r0Mult <= 0 ) then
       ! Invalid value, use default and print warning
       write(preeqObj%io%message, 1000) preeqObj%options%r0Mult, defaultR0Multiplier
       call preeqObj%io%print(2, 3, preeqObj%io%message)
       preeqObj%options%r0Mult = defaultR0Multiplier
    end if



    ! Ensure the number of preeq. particles emitted is within a valid range:
    if ( preeqObj%options%numPreeqType < minNumPreeqType ) then

       ! Must consider for than 'minNumPreeqType' particles for emission
       write(preeqObj%io%message, 2000) preeqObj%options%numPreeqType
       call preeqObj%io%print(2, 3, preeqObj%io%message)
       write(preeqObj%io%message, 2200) minNumPreeqType
       call preeqObj%io%print(2, 3, preeqObj%io%message)

       preeqObj%options%numPreeqType = minNumPreeqType

    else if ( preeqObj%options%numPreeqType > maxNumPreeqType ) then

       ! Cannot consider this many particles for preequilibrium emissio; limit
       write(preeqObj%io%message, 2100) preeqObj%options%numPreeqType
       call preeqObj%io%print(2, 3, preeqObj%io%message)
       write(preeqObj%io%message, 2200) maxNumPreeqType
       call preeqObj%io%print(2, 3, preeqObj%io%message)

       preeqObj%options%numPreeqType = maxNumPreeqType

    end if



    ! Set which level density coefficients to use
    if (   preeqObj%options%levelDenParam /= 11 .and. &
         & preeqObj%options%levelDenParam /= 12 ) then
       write(preeqObj%io%message, 3000) preeqObj%options%levelDenParam
       call preeqObj%io%print(2, 3, preeqObj%io%message)
       write(preeqObj%io%message, 3010) defaultLevelDenParam
       call preeqObj%io%print(2, 3, preeqObj%io%message)
       preeqObj%options%levelDenParam = defaultLevelDenParam
    end if



    ! Set preequilibrium emission probability factor
    if ( preeqObj%options%emissionWidth <= 0 ) then
       write(preeqObj%io%message, 4000) preeqObj%options%emissionWidth
       call preeqObj%io%print(2, 3, preeqObj%io%message)
       write(preeqObj%io%message, 4010) defaultEmissionWidth
       call preeqObj%io%print(2, 3, preeqObj%io%message)
       preeqObj%options%emissionWidth = defaultEmissionWidth
    end if


! ======================================================================
1000 format("Invalid radius parameter (", f5.3, " [f]) detected. Using ", &
          & f5.3, " [fm] instead.")
2000 format("Preequilibrium emission of ", i3, " fragments cannot be considered.")
2100 format("The Preequilibrium class can consider emission of no more than ", &
          & i3, " fragments. This many will be considered.")
2200 format("   Emission of ", i3, " fragments will instead be considered.")
3000 format("Invalid preequilibrium level density flag was used (", i3, ").")
3010 format("   The default level density parameterization flag of ", i3, &
          & " will be used.")
4000 format("An unphysical preequilibrium emission probability width was ", &
          & "cannot be used (", f7.4, ").")
4010 format("   An emission width of ", f7.4, " will instead be used.")
! ======================================================================
  end subroutine validateOptions
