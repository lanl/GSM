
  function coalesDataConstructorMain( kinEnergyGeV, clientOpt ) result(cData)

! ======================================================================
!
! This function is used to construct a results object for clients to use.
! Clients simply pass in the particle bank (using the coalescence particle type)
! and the number of fragments the client wishes to consider for coalescence
!
!
! Written by CMJ, XCP-3, July 2018
! Edited  by CMJ, XCP-3, Jan. 2019
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: real64
    implicit none

    real(real64), intent(in   ) :: kinEnergyGeV
    type(CoalescenceDataOptions), intent(in   ), optional :: clientOpt
    type(coalescenceData)       :: cData

! ======================================================================

    ! Set specific coalescence radii if client wishes
    if ( present(clientOpt) ) then
       cData%options = clientOpt
    end if

    ! Set default values:
    if ( lowerRadiiEBound < kinEnergyGeV .AND. kinEnergyGeV <= upperRadiiEBound ) then
       cData%radii%coalesRadiiDeut  = defaultCRDeut(2)
       cData%radii%coalesRadiiTrit  = defaultCRTrit(2)
       cData%radii%coalesRadiiAlpha = defaultCRAlpha(2)
       cData%radii%coalesRadiiLFrag = defaultCRLFrag(2)
    else
       cData%radii%coalesRadiiDeut  = defaultCRDeut(1)
       cData%radii%coalesRadiiTrit  = defaultCRTrit(1)
       cData%radii%coalesRadiiAlpha = defaultCRAlpha(1)
       cData%radii%coalesRadiiLFrag = defaultCRLFrag(1)
    end if

    ! Set client-specific options for those specified:
    if ( cData%options%coalesRadiiDeut  >= 0 ) cData%radii%coalesRadiiDeut  = cData%options%coalesRadiiDeut
    if ( cData%options%coalesRadiiTrit  >= 0 ) cData%radii%coalesRadiiTrit  = cData%options%coalesRadiiTrit
    if ( cData%options%coalesRadiiAlpha >= 0 ) cData%radii%coalesRadiiAlpha = cData%options%coalesRadiiAlpha
    if ( cData%options%coalesRadiiLFrag >= 0 ) cData%radii%coalesRadiiLFrag = cData%options%coalesRadiiLFrag


    ! Provide values for all squared values:
    cData%radii%coalesRadiiDeutSqrd  = cData%radii%coalesRadiiDeut **2
    cData%radii%coalesRadiiTritSqrd  = cData%radii%coalesRadiiTrit **2
    cData%radii%coalesRadiiAlphaSqrd = cData%radii%coalesRadiiAlpha **2
    cData%radii%coalesRadiiLFragSqrd = cData%radii%coalesRadiiLFrag **2


    ! The results object was constructed successfully
    cData%constructed = .TRUE.

    return
! ======================================================================
  end function coalesDataConstructorMain

! LocalWords:  clientOpt
