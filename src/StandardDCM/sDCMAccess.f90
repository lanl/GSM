
  logical function properlyConstructed(sDCMObj)

! ====================================================================
!
! This function returns to the client a boolean flag that indicates
! if the Standard DCM class was properly constructed
!
!
! Written by CMJ, XCP-3, 01/2019
!
! ====================================================================

    implicit none
    class(StandardDCM), intent(in) :: sDCMObj

! ====================================================================

    properlyConstructed = sDCMObj%constructed

    return
! ====================================================================
  end function properlyConstructed


  function queryOptions( sDCMObj ) result(options)

! ====================================================================
!
! Returns the 'options' object to the client
!
! ====================================================================

    implicit none
    class(StandardDCM), intent(in) :: sDCMObj
    type(sDCMOptions) :: options

! ====================================================================

    options = sDCMObj%options
    return
! ====================================================================
  end function queryOptions


  function queryMolnix( sDCMObj ) result(molObj)

! ====================================================================
!
! Returns the 'molnix' object to the client
!
! ====================================================================

    use molnixClass, only: Molnix

    implicit none
    class(StandardDCM), intent(in) :: sDCMObj
    type(Molnix), pointer :: molObj

! ====================================================================

    molObj => sDCMObj%molnixE
    return
! ====================================================================
  end function queryMolnix


  function queryRNG( sDCMObj ) result(rng)

! ====================================================================
!
! Returns the RNG procedure pointer to the client
!
! ====================================================================

    implicit none
    class(StandardDCM), intent(in) :: sDCMObj
    procedure(RANDOM), pointer :: rng !=> NULL() ! Instantiation not liked by Intel

! ====================================================================

    rng => sDCMObj%rang
    return
! ====================================================================
  end function queryRNG
