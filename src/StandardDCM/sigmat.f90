
  function sigmat8 (sDCM, l, ms, mb, ksi, iks, t)

! ======================================================================
!
!     Choose cross section type and calculate cross section for
!     a given struck-nucleon rest frame kinetic energy.
!     l always 1 for photons; ms always 0; mb = baryon number.
!     mb may be 1 or 2.
!     ksi = 1 for n - n, p - p, pi+ - p & pi- - n.
!     ksi = 2 for n - p, pi+ - n & pi- - p.
!     ksi = 3 for pi0 - n, pi0 - p
!     iks = 0:   total cross section
!     iks = 1:   elastic cross section
!     iks = 2:   pion charge exchange cross section
!     iks = 3:   pion or gamma absorption cross section
!     iks = 4:   neutral pion production cross section
!     iks = 5:   charged pion production cross section
!     iks = 6:   target & projectile isospin change & neutral pion
!                production (9/11/97)
!     iks = 7:   total one-pion production cross section (3/15/99)
!     iks = 8:   delta production (by gamma) cross section
!
!   Called by: CHINEL POINTE TYPINT.
!
!   Calls: QINTS
!   CEM95 written by S. G. Mashnik
!
!   Edited by A. J. Sierk,  LANL  T-2  February-March, 1996.
!   Edited by AJS, July, 1997.
!   Edited by AJS, December, 1998.
!   Corrected by KKG 28.11.01
!   Modified by A. J. Sierk, LANL T-16  October, 2003.
!   Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMParams, only: zro, two
    use standardDCMData,   only: icst, nsicst

    implicit none
    class(StandardDCM), intent(inout) :: sDCM
    integer(int32), intent(in   ) :: l
    integer(int32), intent(in   ) :: ms
    integer(int32), intent(in   ) :: mb
    integer(int32), intent(in   ) :: ksi
    integer(int32), intent(in   ) :: iks
    real(real64),   intent(in   ) :: t
    real(real64)                  :: sigmat8

    integer(int32) ::  ics, js, kns, nsjs

! ======================================================================

!   Factor with which to multiply deuteron absorption cross section to
!   get effective absorption cross section inside nucleus:
!   Corrected pi-0 absorption for factor of 2 historical error.
!   AJS; 12/16/98.
    real(real64), parameter :: absfct = 4.0_real64

!   Factor with which to multiply gamma absorption cross section to
!   get effective absorption cross section inside nucleus:
!   Corrected gamma absorption according to notice at tabulated data.
!   SGM 05/25/03:
    real(real64), parameter :: absfct2 = 5.0_real64

! ======================================================================
!  From bd2:
!     data icst /
!  Translation:total:|elas:   |total: |elas:  |total:  |elas:   |pi+ n or
!              p+p or|p+p or  |       |       |pi- p or|pi- p or|pi- p  |
!              n+n   |n+n     |p+n    |n+p    |pi+ n   |pi+ n   |SCX    |
!    &         210,    211,    220,    221,     120,     121,     122,
!                           | Absorp: |
!  Translation:total:|elas: |pi+ np:pp|p + p->|p + p->|p + n->|p + n->|
!            pi+ p or|pi+ p |pi- pp:np|p p pi0|p n pi+|p + n +|p + p +|
!            pi- n   |pi- n |pi- pn:nn|n + n->|n + n->|pi0    |pi-    |
!                           |pi+ nn:np|n n pi0|p n pi-
!    &         110,    111,    123,    214,     215,     224,    225,
!
!  Transl:   pi+ p->|pi+ p->|pi- p ->|pi- + p->|pi- + p->|gam + p|gam +p|
!         pi0 pi+ p |n +2pi+|2pi0 + n|p + pi- +|n pi- pi+|-> p + |-> n +|
!            pi- n->|pi- n->|pi+ n ->|pi0      |pi+ + n->| pi0   |pi+
!         pi0 pi- n |p +2pi-|2pi0 + p|         |p pi- pi+| (20)   | (21)
!    &         114,    115,    126,     124,      125,    10111,  10112,
!            Absorp:
!  Transl:   gam +  |gam + p |gam + p |gam + p |gam + p |gam + p|gam + p|
!            2N ->  |-> pi+ p|->pi0 p |->n pi+ |-> delta|-> N + |total  |
!   (absorp) 2N (22)|+ pi- 23|+ pi0 24|+ pi0 25|++ + pi-|2 pi 27| (28)  |
!    &       10113,   10115,  10114,   10116    10118,   10117,   10110/

!     data nsicst /
!                   |Absorp:  |  =0! |
!   Key:    pi+ p or|pi+ pn:pp|pi+ p:|pi+ p:|pi+ n:  |pi0 n |pi0 n  |
!            pi- n  |   or    |n 2pi0|N+ 2pi|N + 2pi | or   | or    |
!            SCX=0! |pi- pn:nn|pi- n:|pi- n:|pi- p:  |pi0 p |pi0 p  |
!                             |p 2pi0|N+ 2pi|N + 2pi |total |elastic|
!    &         112,    113,     116,   117,    127,    130,    131,
!                   |Absorp:  |
!   Key:     pi0+p: |pi0+pn:pn|pi0+p:|pi0+p:  |pi0+p:  |pi0+p: |p + p: |
!            pi+ + n|pi0+pp:pp|p 2pi0|p pi+pi-|n pi+pi0|N + 2pi| SCX   |
!            pi0+n: |pi0+np:np|pi0+n:|pi0+n:  |pi0+n:  |pi0+n: |n + n: |
!            pi- + p|pi0+nn:nn|n 2pi0|n pi+pi-|p pi-pi0|N + 2pi| = 0!  |
!    &         132,    133,     134,   135,     136,     137,    212,
!
!   Key:     p + p: |p+p;n+n |p + p->|n + p: |n + p:  |n + p: |n + p->|
!      pi absorption|chg exch|2N + pi| SCX   | Pion   |n p pi0|N N pi |
!            n + n: |+ pi0   |n + n->|       |absorp. |Same as|       |
!             = 0!  |  = 0!  |2N + pi| = 0!  | = 0!   |  224! |       |
!    &         213,    216,     217,   222,     223,     226,    227/

! ======================================================================

!   Start subroutine:
    sigmat8 = zro
    ics = 10000*l + 1000*ms + 100*mb + 10*ksi + iks
    js = 1
    do js = 1,28
       if (ics == icst(js)) then
!   (js will be <= 19, if no photons]).
          sigmat8 = sDCM%qints (t, js)
          if (js == 10) sigmat8 = absfct*sigmat8
          if (js == 22) then
             sigmat8 = absfct2*sigmat8
             if (t < 0.00224d0) sigmat8 = zro
          endif
          if (js == 20 .or. js == 21 .or. js == 28) then
             if (t < 0.15150d0) sigmat8 = zro
          endif
          if (js >= 23 .and. js <= 27) then
             if (t < 0.321d0) sigmat8 = zro
          endif
          return
       endif
    end do
    nsjs = 1
10  if (ics.ne.nsicst(nsjs)) then
       nsjs = nsjs + 1
       if (nsjs < 21) go to 10
    endif
    kns = nsjs
!  1-4: pi+ + p  or  pi- + n  cross section
    if (kns == 1  .or. kns == 3  .or. kns == 14 .or. kns == 15 .or. &
         & kns == 16 .or. kns == 18 .or. kns == 19) then
       sigmat8 = zro
    elseif (kns == 2) then
       sigmat8 = absfct*sDCM%qints (t, 10)
    elseif (kns == 4) then
!  [117] Total pi+ + p or pi- + n pion production cross section:
       sigmat8 = sDCM%qints (t, 15) + sDCM%qints (t, 16)
    elseif (kns == 5) then
!  [127] Total pi- + p or pi+ + n pion production cross section:
       sigmat8 = sDCM%qints (t, 17) + sDCM%qints (t, 18) + sDCM%qints (t, 19)
!  6-13: pi0-p  or  pi0-n  cross sections
    elseif (kns == 6) then
       sigmat8 = (sDCM%qints (t, 8) + sDCM%qints (t, 5))/two
    elseif (kns == 7) then
       sigmat8 = (sDCM%qints (t, 9) + sDCM%qints (t, 6) - sDCM%qints (t, 7))/two
    elseif (kns == 8) then
       sigmat8 = sDCM%qints (t, 7)
    elseif (kns == 9) then
       sigmat8 = 0.5d0*absfct*sDCM%qints (t, 10)
    elseif (kns == 10) then
       sigmat8 = (sDCM%qints (t, 15) + sDCM%qints (t, 18))/two
    elseif (kns == 11) then
       sigmat8 = (sDCM%qints (t, 16) + sDCM%qints (t, 19))/two
    elseif (kns == 12) then
       sigmat8 = sDCM%qints (t, 17)/two
    elseif (kns == 13) then
!  [137] Total pi0+p or pi0+n pion production cross section:
       sigmat8 = (sDCM%qints (t, 15) + sDCM%qints (t, 16) + sDCM%qints (t, 17) + &
            & sDCM%qints (t, 18) + sDCM%qints (t, 19))/two
!  14-17: p-p  or  n-n  cross sections
    elseif (kns == 17) then
!  [217] Total p+p or n+n pion production cross section:
       sigmat8 = sDCM%qints (t, 11) + sDCM%qints (t, 12)
!  18-21: p-n cross sections
    elseif (kns == 20) then
       sigmat8 = sDCM%qints (t, 13)
    elseif (kns == 21) then
!  [227] Total p+n pion-production cross section:
       sigmat8 = sDCM%qints(t, 13) + two*sDCM%qints(t, 14)
    endif
    return

! ======================================================================
  end function sigmat8
