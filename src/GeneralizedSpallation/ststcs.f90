
  subroutine ststcs (gsmObj, a, z, e, ln, bf, ic)

! ======================================================================
!
!   This subroutine keeps statistics on minimum, maximum, and average
!   Atomic number, charge, angular momentum, fission barrier height, and
!   excitation energy of nuclei in various stages of the reaction.
!
!   Written by A. J. Sierk, LANL T-2, January, 1999.
!   (Removed this calculation from PRECOF.)
!   "Last" change: 13-AUG-2003 by NVM
!   Edited by A. J. Sierk, LANL T-16  October, 2003.
!   Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64

    implicit none
    class(GSM),     intent(inout) :: gsmObj
    real(real64),   intent(in   ) :: a
    real(real64),   intent(in   ) :: z
    real(real64),   intent(in   ) :: e
    integer(int32), intent(in   ) :: ln
    real(real64),   intent(in   ) :: bf
    integer(int32), intent(in   ) :: ic

    real(real64) :: el

! ======================================================================

    real(real64) :: eestt, eestsq, aeqtot, aeqsq, zeqtot, zeqsq, &
         & aemin, aemax, eletot, elesq, elemin, elemax, &
         & eestrmn, eestrmx, zemin, zemax
    integer(int32) :: neq
    common /eqidat/  eestt, eestsq, aeqtot, aeqsq, zeqtot, zeqsq, &
         & aemin, aemax, eletot, elesq, elemin, elemax, &
         & eestrmn, eestrmx, zemin, zemax, neq

    real(real64) :: estart, estarsq, estrmn, estrmx, atot, atsq, &
         & ammin, ammax, eltot, elsq, elmmin, elmmax, ztot, &
         & ztsq, zmmin, zmmax, bftot, bfsq, bfmin, bfmax
    common /fisda2/  estart, estarsq, estrmn, estrmx, atot, atsq, &
         & ammin, ammax, eltot, elsq, elmmin, elmmax, ztot, &
         & ztsq, zmmin, zmmax, bftot, bfsq, bfmin, bfmax

    real(real64) :: epstt, epstsq, apqtot, apqsq, zpqtot, zpqsq, &
         & apmin, apmax, elptot, elpsq, elpmin, elpmax, &
         & epstrmn, epstrmx, zpmin, zpmax
    integer(int32) :: npreq
    common /predat/  epstt, epstsq, apqtot, apqsq, zpqtot, zpqsq, &
         & apmin, apmax, elptot, elpsq, elpmin, elpmax, &
         & epstrmn, epstrmx, zpmin, zpmax, npreq

    real(real64) :: erkt, erksq, artot, arsq, zrtot, zrsq, &
         & armin, armax, elrtot, elrsq, elrmin, elrmax, &
         & erkmn, erkmx, zrmin, zrmax
    integer(int32) :: nres
    common /resdat/  erkt, erksq, artot, arsq, zrtot, zrsq, &
         & armin, armax, elrtot, elrsq, elrmin, elrmax, &
         & erkmn, erkmx, zrmin, zrmax, nres

    real(real64) :: arttot, artsq, artmin, artmax, zrttot, zrtsq, & 
         & zrtmin, zrtmax, elrttot, elrtsq, elrtmin, & 
         & elrtmax
    integer(int32) :: nret
    common /retdat/  arttot, artsq, artmin, artmax, zrttot, zrtsq, & 
         & zrtmin, zrtmax, elrttot, elrtsq, elrtmin, & 
         & elrtmax, nret

! ======================================================================

    ! Allow only one thread to enter estimate nuclei statistics at a time
    !$OMP critical

    el = dble(ln)
    if (ic == 1) then
!  Accumulate statistics on average E*, Z, A, L of nucleus after cascade
!  before preequilibrium decay or Fermi break-up.
       npreq = npreq + 1 ! Increment number of pre-equil decays
       apqtot = apqtot + a ! Add to total A
       apqsq = apqsq + a**2 ! sum of all A^2
       apmin = min (apmin, a) ! Smallest A
       apmax = max (apmax, a) ! Largest A
       zpqtot = zpqtot + z 
       zpqsq = zpqsq + z**2
       zpmin = min (zpmin, z)
       zpmax = max (zpmax, z)
       epstt = epstt + e ! Add to total particle energy
       epstsq = epstsq + e**2 ! sum of all E^2
       epstrmn = min (epstrmn, e) ! E min
       epstrmx = max (epstrmx, e) ! E max
       elptot = elptot + el ! Total angular momentum
       elpsq = elpsq + el**2 ! Total angular momentum squared
       elpmin = min (elpmin, el) ! min ang. mom.
       elpmax = max (elpmax, el) ! max ang. mom.
    elseif (ic == 2) then
!   Keep track of nuclei returned with less than 3 MeV of E* after
!   cascade.
       nret = nret + 1
       arttot = arttot + a
       artsq = artsq + a**2
       artmin = min (artmin, a)
       artmax = max (artmax, a)
       zrttot = zrttot + z
       zrtsq = zrtsq + z**2
       zrtmin = min (zrtmin, z)
       zrtmax = max (zrtmax, z)
       elrttot = elrttot + el
       elrtsq = elrtsq + el**2
       elrtmin = min (elrtmin, el)
       elrtmax = max (elrtmax, el)
    elseif (ic == 3) then
!   Nuclei entering the statistical decay phase:
       neq = neq + 1
       aeqtot = aeqtot + a
       aeqsq = aeqsq + a**2
       aemin = min (aemin, a)
       aemax = max (aemax, a)
       zeqtot = zeqtot + z
       zeqsq = zeqsq + z**2
       zemin = min (zemin, z)
       zemax = max (zemax, z)
       eestt = eestt + e
       eestsq = eestsq + e**2
       eestrmn = min (eestrmn, e)
       eestrmx = max (eestrmx, e)
       eletot = eletot + el
       elesq = elesq + el**2
       elemin = min (elemin, el)
       elemax = max (elemax, el)
    elseif (ic == 4) then
!   Residual nuclei:
       nres = nres + 1
       artot = artot + a
       arsq = arsq + a**2
       armin = min (armin, a)
       armax = max (armax, a)
       zrtot = zrtot + z
       zrsq = zrsq + z**2
       zrmin = min (zrmin, z)
       zrmax = max (zrmax, z)
       elrtot = elrtot + el
       elrsq = elrsq + el**2
       elrmin = min (elrmin, el)
       elrmax = max (elrmax, el)
       erkt = erkt + e
       erksq = erksq + e*e
       erkmn = min (erkmn, e)
       erkmx = max (erkmx, e)
    elseif (ic == 5) then
!  Fissioning nuclei:
       estart = estart + e
       estarsq = estarsq + e**2
       estrmn = min(estrmn, e)
       estrmx = max(estrmx, e)
       eltot = eltot + el
       elsq = elsq + el*el
       elmmax = max(el, elmmax)
       elmmin = min(el, elmmin)
       atot = atot + a
       atsq = atsq + a*a
       ammin = min(a, ammin)
       ammax = max(a, ammax)
       ztot = ztot + z
       ztsq = ztsq + z*z
       zmmin = min(z, zmmin)
       zmmax = max(z, zmmax)
       bftot = bftot + bf
       bfsq = bfsq + bf**2
       bfmin = min(bfmin, bf)
       bfmax = max(bfmax, bf)
    else
       write(gsmObj%io%message, 1000) ic
       call gsmObj%io%print(3, 2, gsmObj%io%message)
    endif

    !$OMP end critical

    return
! ======================================================================
1000 format("The 'ststcs' procedure was called with invalid flag (", i4, ").")
! ======================================================================
  end subroutine ststcs
