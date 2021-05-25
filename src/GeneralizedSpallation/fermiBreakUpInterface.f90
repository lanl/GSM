
  subroutine fermiBreakUpInterface (gsmObj, residINC, gsmRxn)

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
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64, int64
    use gsm_params, only: zro, one, thousand, twpi
    use fermiBreakupClass, only : &
         & fermiBreakUpResults, newFermiBreakUpResults, &
         & fermiBreakUpProgeny, fermiBreakUpNucleus

    implicit none
    class(GSM),   intent(inout) :: gsmObj
    class(GSMResidual), intent(in   ) :: residINC
    class(GSMReaction), intent(inout) :: gsmRxn

    integer(int32) :: iafr, ifr, izfr, jfr
    real(real64)   :: cosPhi, cosTheta, phi, totLinMome, &
         & sinPhi, sinTheta, theta, temp

    ! Results object:
    type(fermiBreakUpProgeny), dimension( nint(residINC%numBaryons) ) :: progenyBnk
    type(fermiBreakUpResults) :: fbuResults
    ! Nucleus type:
    type(fermiBreakUpNucleus) :: residual

    ! For photon emission:
    type(GSMResidual) :: photoResidual
    type(GSMResults), pointer :: results => NULL()

! ======================================================================

    ! Set up variables and nucleus:
    residual%numBaryons = residINC%numBaryons
    residual%numProtons = residINC%numProtons
    residual%kinEnergy  = residINC%kinEnergy/thousand
    residual%linearXMom = residINC%linearMom(1)
    residual%linearYMom = residINC%linearMom(2)
    residual%linearZMom = residINC%linearMom(3)
    gsmRxn%outData%ifermi = gsmRxn%outData%ifermi + 1

    ! Create a object for photon emission, if used:
    if ( gsmObj%usePhotonEmission ) then
       photoResidual%numBaryons    = residINC%numBaryons
       photoResidual%numProtons    = residINC%numProtons
       photoResidual%kinEnergy     = residINC%kinEnergy
       photoResidual%linearMom(1)  = residINC%linearMom(1)
       photoResidual%linearMom(2)  = residINC%linearMom(2)
       photoResidual%linearMom(3)  = residINC%linearMom(3)
    end if

    ! Construct the results object:
    fbuResults = newFermiBreakUpResults( progenyBnk )

    ! Perform Fermi Break-Up simulation
    call gsmObj%genModels%fbu%execute (residual, fbuResults)

    ! Now interface results:
    results => gsmRxn%results

    ! Set evaporation/fission flags
    results%info%fusion = zro
    results%info%wf = zro

    ! If fragments created, store data into secondary arrays of program
    progenyLoop: do ifr = 1,fbuResults%numProgeny

       ! Check for available particle storage:
       if ( results%numProgeny > results%maxProgenyM1 ) then
          write(gsmObj%io%message, 1000) fbuResults%numProgeny - ifr + 1
          call gsmObj%io%print(3, 3, gsmObj%io%message)
          exit progenyLoop
       end if

       ! Obtain interim results:
       totLinMome = sqrt( progenyBnk(ifr)%linearXMom**2 + progenyBnk(ifr)%linearYMom**2 + &
            & progenyBnk(ifr)%linearZMom**2 )
       cosTheta = progenyBnk(ifr)%linearZMom/totLinMome
       if ( abs(cosTheta) > one ) cosTheta = sign(one, cosTheta) ! Ensure cos(theta) in limits of [-1, 1]
       sinTheta = sqrt(one - cosTheta**2)
       if (sinTheta > zro) then
          ! Particle has forward movement (get phi component of direction)
          temp = totLinMome*sinTheta
          cosPhi = progenyBnk(ifr)%linearXMom/(temp) ! cos(phi) = px/[sin(theta) ptot]
          sinPhi = progenyBnk(ifr)%linearYMom/(temp) ! sin(phi) = py/[sin(theta) ptot]
       else
          ! Scatter is at 90 degrees (in theta) to incident beam
          cosPhi = one
          sinPhi = zro
       endif
       ! Obtain (phi), ensure within range [0, 2pi]
       phi = atan2 (sinPhi, cosPhi)
       if (phi < zro) phi = twpi + phi
       theta = atan2 (sinTheta, cosTheta)

       ! Obtain GSM fragment label:
       iafr = nint(progenyBnk(ifr)%numBaryons)
       izfr = nint(progenyBnk(ifr)%numProtons)
       jfr = zro
       if (progenyBnk(ifr)%numBaryons <= 4.1d0 .and. progenyBnk(ifr)%numProtons <= 2.1d0) then
          if (iafr == 1 .and. izfr == 0)  jfr = 1 ! n
          if (iafr == 1 .and. izfr == 1)  jfr = 2 ! p
          if (iafr == 2 .and. izfr == 1)  jfr = 3 ! d
          if (iafr == 3 .and. izfr == 1)  jfr = 4 ! t
          if (iafr == 3 .and. izfr == 2)  jfr = 5 ! He-3
          if (iafr == 4 .and. izfr == 2)  jfr = 6 ! He-4
       else
          jfr = thousand*izfr + (iafr - izfr) ! Type for all others
       endif
       if (jfr == 0) jfr = thousand*izfr + (iafr - izfr)

       ! Store progeny in particle bank:
       results%numProgeny = results%numProgeny + 1
       results%progenyBnk(results%numProgeny)%numBaryons = progenyBnk(ifr)%numBaryons
       results%progenyBnk(results%numProgeny)%numProtons = progenyBnk(ifr)%numProtons
       results%progenyBnk(results%numProgeny)%kinEnergy  = progenyBnk(ifr)%kinEnergy / thousand
       results%progenyBnk(results%numProgeny)%restMass   = progenyBnk(ifr)%restMass/thousand
       results%progenyBnk(results%numProgeny)%phi        = phi
       results%progenyBnk(results%numProgeny)%theta      = theta
       results%progenyBnk(results%numProgeny)%sinTheta   = sinTheta
       results%progenyBnk(results%numProgeny)%cosTheta   = cosTheta
       results%progenyBnk(results%numProgeny)%typeID     = jfr
       results%progenyBnk(results%numProgeny)%prodMech   = 1500

       ! Simulate photon emission:
       if ( gsmObj%usePhotonEmission ) then
          call gammaCascade(results%progenyBnk(results%numProgeny), photoResidual)
       end if

    end do progenyLoop

    return
! ======================================================================
1000 format("The GSM progeny array was exceeded. Cannot tally last ", i3, &
          & " fragments.")
! ======================================================================
  end subroutine fermiBreakUpInterface
