
  subroutine  disnmul (gsmObj, fisevent, nn, nnc, nnp, nne, nnprf, nnpof)

! ======================================================================
!
!   This subroutine accumulates the distributions of the multiplicities
!   of neutrons:
!
!     nn   -  total number of neutrons produced in this event;
!     nnc  -  number of neutrons produced in the cascade stage;
!     nnp  -  number of neutrons produced in the preequilibrium stage;
!     nne  -  number of neutrons produced in evaporation (no fission);
!     nnprf-  number of prefission neutrons produced in evaporation;
!     nnpof-  number of postfission neutrons produced in evaporation
!             from fission fragments.
!
!   Written by K. K. Gudima, December, 2004.
!   Edited by A. J. Sierk, LANL T-16, January, 2005.
!   Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use gsm_params, only: one

    implicit none
    class(GSM),     intent(inout) :: gsmObj
    logical,        intent(in   ) :: fisevent
    integer(int32), intent(in   ) :: nn
    integer(int32), intent(in   ) :: nnc
    integer(int32), intent(in   ) :: nnp
    integer(int32), intent(in   ) :: nne
    integer(int32), intent(in   ) :: nnprf
    integer(int32), intent(in   ) :: nnpof

    integer(int32) :: k
    integer(int32), dimension(6) :: n = 0_int32

! ======================================================================

    real(real64) :: disnm
    common /disnmu/ disnm(6,155)

! ======================================================================

    n(1) = nn
    n(2) = nnc
    n(3) = nnp
    n(4) = nne
    n(5) = nnprf
    n(6) = nnpof
    do k=1,6
       if ((k <= 3) .or. (k == 4 .and. .not.fisevent) .or. &
            & (k >= 5 .and. fisevent)) then
          if (n(k) >= 0 .and. n(k) <= 151) then
             disnm(k,n(k)+1) = disnm(k,n(k)+1) + one
             disnm(k,153)    = disnm(k,153) + dble(n(k))
             disnm(k,154)    = disnm(k,154) + dble(n(k)*n(k))
             disnm(k,155)    = disnm(k,155) + one
          endif
       endif
    enddo
    return

! ======================================================================
  end subroutine disnmul
