
  function photoCrossSection (photonEG, e, a0)

! ======================================================================
!
!    "Last" change: 13-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64 

    implicit none
    class(PhotonEventGenerator), intent(inout) :: photonEG
    real(real64),   intent(in   ) :: e
    real(real64),   intent(in   ) :: a0
    real(real64)                  :: photoCrossSection

    real(real64) :: sigkos

! ======================================================================

    if ( abs(a0-one)<divZerLim ) then
       call skosgp (e, sigkos)
    elseif ( abs(a0-two)<divZerLim) then
       call photonEG%skosgd (e, sigkos)
    else
       call photonEG%skosga (e, a0, sigkos)
    endif

    photoCrossSection = sigkos
    return

! ======================================================================
  end function photoCrossSection
