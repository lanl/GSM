
  function properlyConstructed( fbObj ) result(constructed)

! ====================================================================
!
! Returns to client whether or not the FB object was properly constructed
!
! ===================================================================

    implicit none
    class(FissionBarrier), intent(in   ) :: fbObj
    logical                              :: constructed

! ====================================================================

    constructed = fbObj%constructed
    return
! ====================================================================
  end function properlyConstructed


  function queryOptions( fbObj ) result(options)

! ====================================================================
!
! Returns to client the options used by the fission barrier object
!
! ====================================================================

    implicit none
    class(FissionBarrier), intent(in   ) :: fbObj
    type(fissionBarrierOptions)          :: options

! ====================================================================

    options = fbObj%options
    return
! ====================================================================
  end function queryOptions


  function queryMolnix( fbObj ) result(molObj)

! ====================================================================
!
! Returns to client the molnix object used by the fission barrier object
!
! ====================================================================

    use molnixClass, only: Molnix

    implicit none
    class(FissionBarrier), intent(in   ) :: fbObj
    type(Molnix),          pointer       :: molObj

! ====================================================================

    molObj => fbObj%fbMolnix
    return
! ====================================================================
  end function queryMolnix
