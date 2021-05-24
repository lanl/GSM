
  subroutine evapFermiInterface (evapObj, numBaryons, numProtons, &
       & kinEnergyMeV, linearMomx, linearMomy, linearMomz, results)

! ======================================================================
!
!    Fermi break-up calculation of nuclei with A<13 in Preco and Evap
!
!    Called from PRECOF
!
!    Written by K.K. Gudima, 06/23/06
!    Modified by SGM, 07/09/06
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection)
!    Edited by CMJ, XCP-3, 08/2018 (Evap-fermi Interface)
!
! ======================================================================

    use, intrinsic      :: iso_fortran_env, only: int32, real64
    use evaporationParams, only : zro, one, thsn
    use fermiBreakupClass, only : fermiBreakUpProgeny, &
         & fermiBreakUpResults, newFermiBreakUpResults, &
         & fermiBreakUpNucleus


    implicit none
    class(Evaporation), intent(inout) :: evapObj
    real(real64),       intent(in   ) :: numBaryons
    real(real64),       intent(in   ) :: numProtons
    real(real64),       intent(in   ) :: kinEnergyMeV
    real(real64),       intent(in   ) :: linearMomx
    real(real64),       intent(in   ) :: linearMomy
    real(real64),       intent(in   ) :: linearMomz
    type(evaporationResults), intent(inout) :: results

    integer(int32) :: ifr
    real(real64)   :: pfr, pxfr, pyfr, pzfr

    type(fermiBreakUpProgeny), dimension( nint(numBaryons) ) :: fbuProgeny
    type(fermiBreakUpResults) :: fbuResults
    type(fermiBreakUpNucleus) :: fbuNucleus

! ======================================================================

    ! Set up variables
    fbuNucleus%numBaryons = numBaryons
    fbuNucleus%numProtons = numProtons
    fbuNucleus%kinEnergy  = kinEnergyMeV/thsn
    fbuNucleus%linearXMom = linearMomx
    fbuNucleus%linearYMom = linearMomy
    fbuNucleus%linearZMom = linearMomz

    ! Construct FBU results object:
    fbuResults = newFermiBreakUpResults( fbuProgeny )


    ! Perform Fermi Break-Up simulation
    call evapObj%fbuObj%execute (fbuNucleus, fbuResults)


    ! Compound nucleus broke apart; can't fission anymore
    results%fissionProbability = zro

    ! If fragments created, store data into secondary arrays of program
    if (fbuResults%numProgeny > 0) then
       do ifr = 1,fbuResults%numProgeny
          results%numProgeny = results%numProgeny + 1

          ! Ensure within limits
          if ( results%numProgeny >= results%progenyBnkSize ) then
             write(evapObj%io%message,1000)
             call evapObj%io%print(2, 3, evapObj%io%message)
             write(evapObj%io%message,1010) ( results%numProgeny - ifr + 1 )
             call evapObj%io%print(2, 3, evapObj%io%message)
             return
          end if

          ! Obtain normalized momentum of progeny
          pxfr = fbuProgeny(ifr)%linearXMom
          pyfr = fbuProgeny(ifr)%linearYMom
          pzfr = fbuProgeny(ifr)%linearZMom
          pfr  = sqrt( pxfr**2 + pyfr**2 + pzfr**2 )
          if ( pfr < div0Lim ) then
             write(evapObj%io%message, 2000) "95-97"
             call evapObj%io%print(4, 3, evapObj%io%message)
             pfr = div0Lim
          end if

          ! Store into evaporation fragment bank
          results%progenyBnk(results%numProgeny)%numBaryons = fbuProgeny(ifr)%numBaryons
          results%progenyBnk(results%numProgeny)%numProtons = fbuProgeny(ifr)%numProtons
          results%progenyBnk(results%numProgeny)%kinEnergy  = fbuProgeny(ifr)%kinEnergy
          results%progenyBnk(results%numProgeny)%restMass   = &
               & fbuProgeny(ifr)%restMass / thsn
          ! Normalized momentum components
          results%progenyBnk(results%numProgeny)%linearMomX = ( pxfr / pfr )
          results%progenyBnk(results%numProgeny)%linearMomY = ( pyfr / pfr )
          results%progenyBnk(results%numProgeny)%linearMomZ = ( pzfr / pfr )
          results%progenyBnk(results%numProgeny)%physicsFlag = fermiBreakUpFlag

          if ( evapObj%usePhotoEmission ) then
             call evapObj%photonEmission( results%progenyBnk( results%numProgeny ), &
                  &  numBaryons, numProtons, kinEnergyMeV )
          end if


       enddo
    endif

    return

! ======================================================================
1000 format("The evaporation fragment bank was exceeded in the ", &
          & "Fermi Break-up routine.")
1010 format("   The ", i3, " remaining fragments will not be ", &
          & "tallied for this event.")
2000 format("Divide by zero error prevented in ", &
          & "'fermiEvapInterface.f90', line ", A)
! ======================================================================

  end subroutine evapFermiInterface
