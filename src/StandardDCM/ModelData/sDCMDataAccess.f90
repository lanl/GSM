
  function checkIndex( dataObj, i ) result(closestBound)

! ==============================================================================
!
! This function returns an index for a valid entry in the target object arrays.
! These arrays include radii and normalized radii, nucleon densities, and nucleon
! Fermi momenta.
!
!
! Written by CMJ, XCP-3, 01/2019
!
! ==============================================================================

    use, intrinsic:: iso_fortran_env, only: int32

    implicit none
    class(StandardDCMData), intent(inout) :: dataObj
    integer(int32), intent(in   ) :: i
    integer(int32)                :: closestBound

! ==============================================================================

    ! Assume closest bound is what was passed in
    closestBound = i

    ! Check if the element can exist, if not, move to the closest existing bound
    if ( i < minZonesAllowed ) then
       write(dataObj%io%message, 1000) i, minZonesAllowed
       call dataObj%io%print(0, 3, dataObj%io%message)
       closestBound = minZonesAllowed
    else if ( i > maxZonesAllowed ) then
       write(dataObj%io%message, 1000) i, maxZonesAllowed
       call dataObj%io%print(0, 3, dataObj%io%message)
       closestBound = maxZonesAllowed
    end if

    return
! ==============================================================================
1000 format("Standard DCM data array exceeded - out of bounds (element ", i3, &
          & "). Using closest valid index (", i3, ").")
! ==============================================================================
  end function checkIndex



  function zoneBoundR( dataObj, i ) result(boundR)

! ==============================================================================
!
! This function returns to the client the value of the 'zoneBoundR' variable.
!
!
! Written by CMJ, XCP-3, 01/2019
!
! ==============================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    class(StandardDCMData), intent(inout) :: dataObj
    integer(int32),         intent(in   ) :: i

    integer(int32) :: indx
    real(real64)   :: boundR

! ==============================================================================

    indx = dataObj%checkIndex( i )

    boundR = dataObj%target%zoneBoundR( indx )

    return
! ==============================================================================
  end function zoneBoundR



  function zoneBRFrac( dataObj, i ) result(boundRF)

! ==============================================================================
!
! This function returns to the client the value of the 'zoneBRFrac' variable.
!
!
! Written by CMJ, XCP-3, 01/2019
!
! ==============================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    class(StandardDCMData), intent(inout) :: dataObj
    integer(int32),         intent(in   ) :: i

    integer(int32) :: indx
    real(real64)   :: boundRF

! ==============================================================================

    indx = dataObj%checkIndex( i )

    boundRF = dataObj%target%zoneBRFrac( indx )

    return
! ==============================================================================
  end function zoneBRFrac



  function protonDensity( dataObj, i ) result(density)

! ==============================================================================
!
! This function returns to the client the value of the 'protonDensity' variable.
!
!
! Written by CMJ, XCP-3, 01/2019
!
! ==============================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    class(StandardDCMData), intent(inout) :: dataObj
    integer(int32),         intent(in   ) :: i

    integer(int32) :: indx
    real(real64)   :: density

! ==============================================================================

    indx = dataObj%checkIndex( i )

    density = dataObj%target%protonDensity( indx )

    return
! ==============================================================================
  end function protonDensity



  function neutronDensity( dataObj, i ) result(density)

! ==============================================================================
!
! This function returns to the client the value of the 'neutronDensity' variable.
!
!
! Written by CMJ, XCP-3, 01/2019
!
! ==============================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    class(StandardDCMData), intent(inout) :: dataObj
    integer(int32),         intent(in   ) :: i

    integer(int32) :: indx
    real(real64)   :: density

! ==============================================================================

    indx = dataObj%checkIndex( i )

    density = dataObj%target%neutronDensity( indx )

    return
! ==============================================================================
  end function neutronDensity



  function coulombPote( dataObj, i ) result(coulV)

! ==============================================================================
!
! This function returns to the client the value of the 'coulombPote' variable.
!
!
! Written by CMJ, XCP-3, 01/2019
!
! ==============================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    class(StandardDCMData), intent(inout) :: dataObj
    integer(int32),         intent(in   ) :: i

    integer(int32) :: indx
    real(real64)   :: coulV

! ==============================================================================

    indx = dataObj%checkIndex( i )

    coulV = dataObj%target%coulombPote( indx )

    return
! ==============================================================================
  end function coulombPote



  function protFermiMom( dataObj, i ) result(fProtMom)

! ==============================================================================
!
! This function returns to the client the value of the 'protFermiMom' variable.
!
!
! Written by CMJ, XCP-3, 01/2019
!
! ==============================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    class(StandardDCMData), intent(inout) :: dataObj
    integer(int32),         intent(in   ) :: i

    integer(int32) :: indx
    real(real64)   :: fProtMom

! ==============================================================================

    indx = dataObj%checkIndex( i )

    fProtMom = dataObj%target%protFermiMom( indx )

    return
! ==============================================================================
  end function protFermiMom



  function neutFermiMom( dataObj, i ) result(fNeutMom)

! ==============================================================================
!
! This function returns to the client the value of the 'neutFermiMom' variable.
!
!
! Written by CMJ, XCP-3, 01/2019
!
! ==============================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    class(StandardDCMData), intent(inout) :: dataObj
    integer(int32),         intent(in   ) :: i

    integer(int32) :: indx
    real(real64)   :: fNeutMom

! ==============================================================================

    indx = dataObj%checkIndex( i )

    fNeutMom = dataObj%target%neutFermiMom( indx )

    return
! ==============================================================================
  end function neutFermiMom



  function numBaryons ( dataObj ) result(aNumber)

! ==============================================================================
!
! This function returns the number of baryons contained by the target that the
! Standard DCM Data object is established for.
!
!
! Written by CMJ, XCP-3, 01/2019
!
! ==============================================================================

    use, intrinsic:: iso_fortran_env, only: real64

    implicit none
    class(StandardDCMData), intent(in   ) :: dataObj
    real(real64) :: aNumber

! ==============================================================================

    aNumber = dataObj%target%numBaryons
    return
! ==============================================================================
  end function numBaryons


  function numProtons ( dataObj ) result(zNumber)

! ==============================================================================
!
! This function returns the number of protons contained by the target that the
! Standard DCM Data object is established for.
!
!
! Written by CMJ, XCP-3, 01/2019
!
! ==============================================================================

    use, intrinsic:: iso_fortran_env, only: real64

    implicit none
    class(StandardDCMData), intent(in   ) :: dataObj
    real(real64) :: zNumber

! ==============================================================================

    zNumber = dataObj%target%numProtons
    return
! ==============================================================================
  end function numProtons


  function numZones ( dataObj ) result(nZones)

! ==============================================================================
!
! This function returns the number of zones contained by the target that the
! Standard DCM Data object is established for.
!
!
! Written by CMJ, XCP-3, 01/2019
!
! ==============================================================================

    use, intrinsic:: iso_fortran_env, only: int32

    implicit none
    class(StandardDCMData), intent(in   ) :: dataObj
    integer(int32) :: nZones

! ==============================================================================

    nZones = dataObj%options%numZones
    return
! ==============================================================================
  end function numZones


  function aTargThrd ( dataObj ) result(athrd)

! ==============================================================================
!
! This function returns the aTarg**(1/3) for the Standard DCM Data object.
!
!
! Written by CMJ, XCP-3, 01/2019
!
! ==============================================================================

    use, intrinsic:: iso_fortran_env, only: real64

    implicit none
    class(StandardDCMData), intent(in   ) :: dataObj
    real(real64)                          :: athrd

! ==============================================================================

    athrd = dataObj%target%aTargThrd
    return
! ==============================================================================
  end function aTargThrd


  function geomCrossSection ( dataObj ) result(geomSigma)

! ==============================================================================
!
! This function returns the geometric cross section of the target that the
! Standard DCM Data object is established for.
!
!
! Written by CMJ, XCP-3, 01/2019
!
! ==============================================================================

    use, intrinsic:: iso_fortran_env, only: real64

    implicit none
    class(StandardDCMData), intent(in   ) :: dataObj
    real(real64)                          :: geomSigma

! ==============================================================================

    geomSigma = dataObj%target%geomCrossSection
    return
! ==============================================================================
  end function geomCrossSection


  function pionPote ( dataObj ) result(pionPotential)

! ==============================================================================
!
! This function returns the geometric cross section of the target that the
! Standard DCM Data object is established for.
!
!
! Written by CMJ, XCP-3, 01/2019
!
! ==============================================================================

    use, intrinsic:: iso_fortran_env, only: real64

    implicit none
    class(StandardDCMData), intent(in   ) :: dataObj
    real(real64)                          :: pionPotential

! ==============================================================================

    pionPotential = dataObj%target%pionPote
    return
! ==============================================================================
  end function pionPote


  function getSepEnergy ( dataObj ) result(sepEnergy)

! ==============================================================================
!
! This function returns the geometric cross section of the target that the
! Standard DCM Data object is established for.
!
!
! Written by CMJ, XCP-3, 01/2019
!
! ==============================================================================

    use, intrinsic:: iso_fortran_env, only: real64

    implicit none
    class(StandardDCMData), intent(in   ) :: dataObj
    real(real64)                          :: sepEnergy

! ==============================================================================

    sepEnergy = dataObj%target%sepEnergy
    return
! ==============================================================================
  end function getSepEnergy


  subroutine setSepEnergy( dataObj, newSepEnergy )

! ==============================================================================
!
! This function returns the geometric cross section of the target that the
! Standard DCM Data object is established for.
!
!
! Written by CMJ, XCP-3, 01/2019
!
! ==============================================================================

    use, intrinsic:: iso_fortran_env, only: real64

    implicit none
    class(StandardDCMData), intent(inout) :: dataObj
    real(real64),           intent(in   ) :: newSepEnergy

! ==============================================================================

    dataObj%target%sepEnergy = newSepEnergy
    return
! ==============================================================================
  end subroutine setSepEnergy


  function properlyConstructedObject(dataObj) result(constructed)

! ==============================================================================
!
! Returns flag stating if the data object is constructed or not
!
! ==============================================================================

    implicit none
    class(StandardDCMData), intent(in   ) :: dataObj
    logical :: constructed

! ==============================================================================

    constructed = dataObj%objectConstructed

    return
! ==============================================================================
  end function properlyConstructedObject


  function queryOptions( dataObj ) result(options)

! ==============================================================================
!
! Returns the options contained by the sDCM Data object
!
! ==============================================================================

    implicit none
    class(StandardDCMData), intent(in   ) :: dataObj
    type(sDCMDataOptions) :: options

! ==============================================================================

    options = dataObj%options
    return
! ==============================================================================
  end function queryOptions
