
  subroutine setPhotonEmission( clientPhotonEmission )

! ====================================================================
!
! Allows clients to utilize a SINGLE photon emission procedure for the
! entirety of GSM
!
! NOTE: If photon emission is NOT specified prior to GSM object
!       construction, it will NOT be used.
!
!
! Written by CMJ, XCP-3 (04/2019)
!
! ====================================================================

    implicit none
    procedure(PHOTOEMISSION), intent(in), pointer :: clientPhotonEmission

    ! To handle messagging temporarily:
    type(GSMIO) :: io

! ====================================================================

    if ( associated(clientPhotonEmission) ) then
       useGammaCascade = .TRUE.
       gammaCascade => clientPhotonEmission
    else
       write(io%message, 1000)
       call io%print(3, 2, io%message)
    end if

    return
! ====================================================================
1000 format("The GSM-provided photon emission procedure is not ", &
          & "associated and will not be used.")
! ====================================================================
  end subroutine setPhotonEmission


  subroutine preeqPhotonEmission(preeqFrag, residual)

! ====================================================================
!
! Interfaces from preequilibrium emission procedure to a GSM-level
! one for client access to photon emission
!
! ====================================================================

    use gsm_params, only: thousand
    use preequilibriumClass, only: preequilibriumFragment, residualNucleus

    implicit none
    type(preequilibriumFragment), intent(in) :: preeqFrag
    type(residualNucleus),        intent(in) :: residual

    type(GSMProgeny)  :: progeny
    type(GSMResidual) :: parentNucleus

! ====================================================================

    ! Obtain progeny:
    progeny%numBaryons = preeqFrag%numBaryons
    progeny%numProtons = preeqFrag%numProtons
    progeny%kinEnergy  = preeqFrag%kinEnergy / thousand
    progeny%restMass   = preeqFrag%restMass / thousand
    progeny%phi        = preeqFrag%phi
    progeny%theta      = preeqFrag%theta
    progeny%sinTheta   = sin( progeny%theta )
    progeny%cosTheta   = cos( progeny%theta )
    progeny%typeID     = 0
    progeny%prodMech   = 100

    ! Obtain parent nucleus:
    ! NOTE: All units match
    parentNucleus%numBaryons = residual%numBaryons
    parentNucleus%numProtons = residual%numProtons
    parentNucleus%kinEnergy  = residual%kinEnergy
    parentNucleus%linearMom(1) = residual%linearMomX
    parentNucleus%linearMom(2) = residual%linearMomY
    parentNucleus%linearMom(3) = residual%linearMomZ
    parentNucleus%angularMom(1) = residual%angularMom(1)
    parentNucleus%angularMom(2) = residual%angularMom(2)
    parentNucleus%angularMom(3) = residual%angularMom(3)

    ! Use photon emission procedure:
    ! NOTE: This is a redundant check (see preeq. model construction)!
    if ( useGammaCascade ) then
       call gammaCascade(progeny, parentNucleus)
    end if

    return
! ====================================================================
  end subroutine preeqPhotonEmission


  subroutine evapPhotonEmission(evapFrag, residA, residZ, residKE)

! ====================================================================
!
! Interfaces from preequilibrium emission procedure to a GSM-level
! one for client access to photon emission
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use gsm_params, only: thousand
    use evaporationClass, only: evaporationFragment

    implicit none
    type(evaporationFragment), intent(in) :: evapFrag
    real(real64),              intent(in) :: residA
    real(real64),              intent(in) :: residZ
    real(real64),              intent(in) :: residKE

    type(GSMProgeny)  :: progeny
    type(GSMResidual) :: parentNucleus

! ====================================================================

    ! Obtain progeny:
    progeny%numBaryons = evapFrag%numBaryons
    progeny%numProtons = evapFrag%numProtons
    progeny%kinEnergy  = evapFrag%kinEnergy / thousand
    progeny%restMass   = evapFrag%restMass
    progeny%phi        = atan( evapFrag%linearMomY / evapFrag%linearMomX )
    progeny%theta      = acos( evapFrag%linearMomZ / 1.0_real64 )
    progeny%sinTheta   = sin( progeny%theta )
    progeny%cosTheta   = cos( progeny%theta )
    progeny%typeID     = 0
    progeny%prodMech   = 1000

    ! Obtain parent nucleus:
    ! NOTE: All units match
    parentNucleus%numBaryons = residA
    parentNucleus%numProtons = residZ
    parentNucleus%kinEnergy  = residKE

    ! Use photon emission procedure:
    ! NOTE: This is a redundant check (see evap. model construction)!
    if ( useGammaCascade ) then
       call gammaCascade(progeny, parentNucleus)
    end if

    return
! ====================================================================
  end subroutine evapPhotonEmission
