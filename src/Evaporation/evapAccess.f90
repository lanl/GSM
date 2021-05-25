
! ====================================================================
!
! This file contains all procedures clients may use to query the
! evaporation object.
!
!
! Written by CMJ, XCP03 (03/2019)
!
! ====================================================================

  function properlyConstructed(evapObj) result(constructed)

! ====================================================================
!
! Returns whether or not the class is properly constructed
!
! ====================================================================

    implicit none
    class(Evaporation), intent(in   ) :: evapObj
    logical :: constructed

! ====================================================================

    constructed = evapObj%constructed

    return
! ====================================================================
  end function properlyConstructed


  function queryOptions(evapObj) result(options)

! ====================================================================
!
! Returns the objects options set
!
! ====================================================================

    implicit none
    class(Evaporation), intent(in   ) :: evapObj
    type(evaporationOptions) :: options

! ====================================================================

    options = evapObj%options

    return
! ====================================================================
  end function queryOptions


  function queryRNG(evapObj) result(rngPtr)

! ====================================================================
!
! Returns the object's RNG
!
! ====================================================================

    implicit none
    class(Evaporation), intent(in   ) :: evapObj
    procedure(RANDOM), pointer :: rngPtr

! ====================================================================

    rngPtr => evapObj%rang

    return
! ====================================================================
  end function queryRNG


  function queryMolnix(evapObj) result(molObj)

! ====================================================================
!
! Returns the object's Molnix object
!
! ====================================================================

    use molnixClass, only: Molnix

    implicit none
    class(Evaporation), intent(in   ) :: evapObj
    class(Molnix), pointer :: molObj

! ====================================================================

    molObj => evapObj%evapMolnix

    return
! ====================================================================
  end function queryMolnix


  function queryData(evapObj) result(evapData)

! ====================================================================
!
! Returns the object's evaporation data object
!
! ====================================================================

    use evaporationDataClass, only: EvaporationData

    implicit none
    class(Evaporation), intent(in   ) :: evapObj
    class(EvaporationData), pointer :: evapData

! ====================================================================

    evapData => evapObj%data

    return
! ====================================================================
  end function queryData


  function queryFermiBreakUp(evapObj) result(fbuObj)

! ====================================================================
!
! Returns the object's Fermi Breakup object
!
! ====================================================================

    use fermiBreakUpClass, only: FermiBreakup

    implicit none
    class(Evaporation), intent(in   ) :: evapObj
    type(FermiBreakup) :: fbuObj

! ====================================================================

    fbuObj = evapObj%fbuObj

    return
! ====================================================================
  end function queryFermiBreakUp


  function queryPhotoEmissionUse(evapObj) result(useEmission)

! ====================================================================
!
! Returns a flag indicating if photon emission is used or not
!
! ====================================================================

    implicit none
    class(Evaporation), intent(in   ) :: evapObj
    logical :: useEmission

! ====================================================================

    useEmission = evapObj%usePhotoEmission

    return
! ====================================================================
  end function queryPhotoEmissionUse


  function queryPhotoEmission(evapObj) result(gammaEmission)

! ====================================================================
!
! Returns the object's photon emission procedure pointer
!
! ====================================================================

    implicit none
    class(Evaporation), intent(in   ) :: evapObj
    procedure(PHOTOEMISSION), pointer :: gammaEmission

! ====================================================================

    gammaEmission => evapObj%photonEmission

    return
! ====================================================================
  end function queryPhotoEmission
