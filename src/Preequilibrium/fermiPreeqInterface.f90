
  subroutine fermiPreeqInterface (preeqObj, residual, results)

! ======================================================================
!
!    Fermi break-up calculation of nuclei with A<13 in Preco and Evap
!
!    Called from EQUILIBRATE
!
!    Written by K.K. Gudima, 06/23/06
!    Modified by SGM, 07/09/06
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection)
!    Edited by CMJ, XCP-3, July 2018 (included in preeq. class)
!
! ======================================================================

    use, intrinsic      :: iso_fortran_env, only: int32, real64
    use preequilibriumParams,   only: zro, one, thousand, twpi
    use fermiBreakupClass,      only: &
         & fermiBreakUpProgeny, &
         & fermiBreakUpresults, newFermiBreakUpResults, &
         & fermiBreakUpNucleus


    implicit none
    class(Preequilibrium), intent(inout) :: preeqObj
    type(residualNucleus), intent(in   ) :: residual
    type(preequilibriumResults), intent(inout) :: results

    integer(int32) :: ifr
    real(real64)   :: cffr, ctfr, fifr, pfr, pxfr, pyfr, pzfr, sffr, &
         & stfr, tefr, temp

    type(fermiBreakUpProgeny), dimension( nint(residual%numBaryons) ) :: fbuProgeny
    type(fermiBreakUpResults) :: fbuResults
    type(fermiBreakUpNucleus) :: fbuNucleus

! ======================================================================

    ! Set up variables
    fbuNucleus%numBaryons = residual%numBaryons
    fbuNucleus%numProtons = residual%numProtons
    fbuNucleus%kinEnergy  = residual%thermEnergy/thousand
    fbuNucleus%linearXMom = residual%linearMomX
    fbuNucleus%linearYMom = residual%linearMomY
    fbuNucleus%linearZMom = residual%linearMomZ


    ! Construct the FBU Results:
    fbuResults = newFermiBreakUpResults( fbuProgeny )


    ! Perform Fermi Break-Up simulation
    call preeqObj%fbuObj%execute (fbuNucleus, fbuResults)


    ! If fragments created, store data into secondary arrays of program
    if (fbuResults%numProgeny > 0) then
       do ifr = 1, fbuResults%numProgeny

          ! Obtain momentums
          pxfr = fbuProgeny(ifr)%linearXMom
          pyfr = fbuProgeny(ifr)%linearYMom
          pzfr = fbuProgeny(ifr)%linearZMom
          pfr  = sqrt(pxfr**2 + pyfr**2 + pzfr**2) ! Total momentum
          if (pfr < div0Lim .and. pfr > -div0Lim) then
             pfr = div0Lim
             write(preeqObj%io%message,1000) "69"
             call preeqObj%io%print(4, 3, preeqObj%io%message)
          end if


          ! Obtain theta, phi (based on momentum vector)
          ctfr = pzfr/pfr ! cos(theta) = pz/ptot (for this fragment)
          stfr = sqrt(abs(one - ctfr**2)) ! sin(theta) = Sqrt[ 1 - cos(theta)^2]
          if (stfr > zro) then
             ! Particle has forward movement (get phi component of direction)
             temp = pfr*stfr
             if (temp < div0Lim .and. temp > -div0Lim) then
                temp = div0Lim
                write(preeqObj%io%message,1000) "76, 77"
                call preeqObj%io%print(4, 3, preeqObj%io%message)
             end if
             cffr = pxfr/(temp) ! cos(phi) = px/[sin(theta) ptot]
             sffr = pyfr/(temp) ! sin(phi) = py/[sin(theta) ptot]
          else
             ! Scatter is at 90 degrees (in theta) to incident beam
             cffr = one
             sffr = zro
          endif


          ! Obtain (phi), ensure within range [0, 2pi]
          fifr = atan2 (sffr, cffr)
          if (fifr < zro) fifr = twpi + fifr
          ! Obtain (theta)
          tefr = atan2 (stfr, ctfr)


          ! Add fragment to list of secondaries
          results%numProgeny = results%numProgeny + 1


          ! Ensure within limits
          if ( results%numProgeny > results%maxProgeny ) then
             ! Preeq. array will be exceeded; print warning and return
             write(preeqObj%io%message,2000)
             call preeqObj%io%print(1, 3, preeqObj%io%message)
             write(preeqObj%io%message,2100) &
                  & results%numProgeny - results%maxProgeny
             call preeqObj%io%print(1, 3, preeqObj%io%message)

             ! Remove progeny from array (unable to store)
             results%numProgeny = results%numProgeny - 1

             return
          end if


          results%progenyBnk(results%numProgeny)%numBaryons = fbuProgeny(ifr)%numBaryons
          results%progenyBnk(results%numProgeny)%numProtons = fbuProgeny(ifr)%numProtons
          results%progenyBnk(results%numProgeny)%kinEnergy  = fbuProgeny(ifr)%kinEnergy
          results%progenyBnk(results%numProgeny)%restMass   = fbuProgeny(ifr)%restMass
          results%progenyBnk(results%numProgeny)%theta      = tefr
          results%progenyBnk(results%numProgeny)%phi        = fifr
          results%progenyBnk(results%numProgeny)%origin     = fermiBreakUpProgenyType
          

          ! Check for photon emission:
          if ( preeqObj%usePhotonEmission ) then
             call preeqObj%photonEmission( results%progenyBnk( results%numProgeny ), &
                  & residual )
          end if

       enddo
    endif

    return

! ======================================================================
1000 format("Divide by zero error prevented in ", &
          & "'fermiInterface.f90', line ", A)
2000 format("The Preequilibrium Fragment array was exceeded during the ", &
          & "Fermi Break-up process.")
2100 format("   The ", i3, " remaining fragments will not be ", &
          & "tallied for this event.")
! ======================================================================
  end subroutine fermiPreeqInterface
