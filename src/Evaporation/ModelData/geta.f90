
  function geta (alev, e, iz, in, is)

! ======================================================================
!
! ***************This routine was originally in the LAHET code**********
!  GETA
!      Calculate level density parameter
! ======================================================================
! <variables>
!     iz   :   charge  of nucleus                                 (IN)
!     in   :   neutron number of nucleus                          (IN)
!     e    :   excitation energy of nucleus                       (IN)
! ======================================================================
!
!    "Last" change: 13-AUG-2003 by NVMokhov
!    Modified by A. J. Sierk, LANL T-16, September, 2003.
!    Modified by K. K. Gudima, November, 2004.
!    Edited by A. J. Sierk, January, 2005.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use evaporationParams,      only: zro, one, eighth
    use evaporationFissionData, only: inn, iiz, aMax, isn, isz, &
         & amean, con, paire0, st0, &
         & levelDensityGCCIFlag, levelDensitySpecFlag

    implicit none
    real(real64),   intent(in   ) :: alev
    real(real64),   intent(in   ) :: e
    integer(int32), intent(in   ) :: iz
    integer(int32), intent(in   ) :: in
    integer(int32), intent(  out) :: is
    real(real64)                  :: geta

    integer(int32) :: ia, ia0
    real(real64)   :: aa, deldel, fa, fu, gu, u

! ======================================================================

    real(real64), parameter, dimension(aMax) :: fa0 = &
       & [   ( (eighth),               ia0 =   1,  25), &
       &     ( (amean(ia0)/dble(ia0)), ia0 =  26, 240), &
       &     ( ((amean(240) + (2.5d0 - 0.1d0*amean(240))* &
       &        (dble(ia0) - 240.d0))/dble(ia0) &
       &                        ),     ia0 = 241, aMax)   ]

! ======================================================================

    ia = in + iz
    aa = dble(ia)
    if (alev /= levelDensityGCCIFlag) then
!   a = A/a0 ...simple level density parameter:
       if ( alev == levelDensitySpecFlag ) then
          geta = aa*eighth
       else
          geta = aa/alev
       endif
       return
    endif

! Gilbert-Cameron-Cook-Ignatyuk (GCCI) level density parameter:
    if (iz > iiz .or. in > inn) then
       geta = aa*eighth
       return
    endif
    if (in >= 9 .and. iz >= 9) then
       if (isz(iz).or.isn(in)) then
          is = 2
       else
          is = 1
       endif
       fa = (9.17d-3*st0(iz,in) + con(is))
    else
       fa = fa0(ia)
    endif

    deldel = paire0(iz, in)
    u = max(5.d-2*(e - deldel), 1.d-5)
    fu = exp(-u)
    gu = (one - fu)/u
    geta = aa*(fa*gu + (one - gu)*(0.1375d0 - 8.36d-5*aa))

    return
! ======================================================================
  end function geta
