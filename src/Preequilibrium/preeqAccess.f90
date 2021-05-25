
! ====================================================================
!
! This file contains all methods by which clients can query the
! preequilibrium object.
!
!
! Written by CMJ, XCP-3 (03/2019)
!
! ====================================================================

  function properlyConstructed(preeqObj) result(constructed)

! ====================================================================
!
! The procedure returns the construction state of the class
!
! ====================================================================

    implicit none
    class(Preequilibrium), intent(in   ) :: preeqObj
    logical :: constructed

! ====================================================================

    constructed = preeqObj%constructed

    return
! ====================================================================
  end function properlyConstructed


  function queryOptions(preeqObj) result(options)

! ====================================================================
!
! The procedure returns the options of the class
!
! ====================================================================

    implicit none
    class(Preequilibrium), intent(in   ) :: preeqObj
    type(preequilibriumOptions) :: options

! ====================================================================

    options = preeqObj%options

    return
! ====================================================================
  end function queryOptions


  function queryMolnix(preeqObj) result(molObj)

! ====================================================================
!
! The procedure returns the Molnix pointer used by the object
!
! ====================================================================

    use molnixClass, only: Molnix

    implicit none
    class(Preequilibrium), intent(in   ) :: preeqObj
    type(Molnix), pointer :: molObj

! ====================================================================

    molObj => preeqObj%molEnergy

    return
! ====================================================================
  end function queryMolnix


  function queryFissionBarrier(preeqObj) result(fissBObj)

! ====================================================================
!
! The procedure returns the fission barrier object pointer used by the class
!
! ====================================================================

    use fissionBarrierClass, only: FissionBarrier

    implicit none
    class(Preequilibrium), intent(in   ) :: preeqObj
    type(FissionBarrier), pointer :: fissBObj

! ====================================================================

    fissBObj => preeqObj%fissBarr

    return
! ====================================================================
  end function queryFissionBarrier


  function queryRNG(preeqObj) result(rang)

! ====================================================================
!
! The procedure returns the RNG procedure pointer used by the class
!
! ====================================================================

    implicit none
    class(Preequilibrium), intent(in   ) :: preeqObj
    procedure(RANDOM), pointer :: rang

! ====================================================================

    rang => preeqObj%rng

    return
! ====================================================================
  end function queryRNG


  function queryFermiBreakUp(preeqObj) result(fbuObj)

! ====================================================================
!
! The procedure returns the FermiBreakup model used by the class
!
! ====================================================================

    use fermiBreakUpClass, only: FermiBreakup

    implicit none
    class(Preequilibrium), intent(in   ) :: preeqObj
    type(FermiBreakup) :: fbuObj

! ====================================================================

    fbuObj = preeqObj%fbuObj

    return
! ====================================================================
  end function queryFermiBreakUp


  function queryPhotonEmissionUse(preeqObj) result(useEmission)

! ====================================================================
!
! The procedure returns the flag indicated use of photon emission used by the class
!
! ====================================================================

    implicit none
    class(Preequilibrium), intent(in   ) :: preeqObj
    logical :: useEmission


! ====================================================================

    useEmission = preeqObj%usePhotonEmission

    return
! ====================================================================
  end function queryPhotonEmissionUse


  function queryPhotonEmission(preeqObj) result(emissionProcedure)

! ====================================================================
!
! The procedure returns the photon emission procedure pointer used by the class
!
! ====================================================================

    implicit none
    class(Preequilibrium), intent(in   ) :: preeqObj
    procedure(PHOTOEMISSION), pointer :: emissionProcedure

! ====================================================================

    emissionProcedure => preeqObj%photonEmission

    return
! ====================================================================
  end function queryPhotonEmission


  function queryData(preeqObj) result(dataObj)

! ====================================================================
!
! The procedure returns the  used by the class
!
! ====================================================================

    use preequilibriumDataClass, only: PreequilibriumData

    implicit none
    class(Preequilibrium), intent(in   ) :: preeqObj
    type(PreequilibriumData), pointer :: dataObj

! ====================================================================

    dataObj => preeqObj%preeqData

    return
! ====================================================================
  end function queryData
