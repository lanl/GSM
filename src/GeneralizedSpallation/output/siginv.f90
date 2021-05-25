
  function sinvla (gsmObj, j, a, z, e)

! ======================================================================
!
!     Written in 2001 by K.K. Gudima
!     Edited in 2001 by S.G. Mashnik
!     Modified in October, 2003 by A. J. Sierk
!     Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================
!
!  Calculates the inelastic cross section (mb) for interaction of the
!  particle "j" of energy e (MeV) with the nucleus (a,z) using NASA
!  systematics (NIM B 117 (1996) 347) for all charged particles and
!  neutrons with energy above the peak in the cross section, and
!  Kalbach systematics (J. Phys. G: Nucl. Phys., 24 (1998) 847) for
!  lower energy neutrons;
!
!  j = 1, 2, 3, 4, 5, 6 means n, p, d, t, He3, and He4
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use gsm_params, only: zro, one
    use generalizedSpallationData, only: enmax, scalf

    implicit none
    class(GSM),     intent(inout) :: gsmObj
    integer(int32), intent(in   ) :: j
    real(real64),   intent(in   ) :: a
    real(real64),   intent(in   ) :: z
    real(real64),   intent(in   ) :: e
    real(real64)                  :: sinvla

    integer(int32) :: iat
    real(real64)   :: ea, enm, scf

! ======================================================================

    real(real64), dimension( 6 ), parameter :: &
         & ap = [ 1.0_real64, 1.0_real64, 2.0_real64, 3.0_real64, 3.0_real64, 4.0_real64 ], &
         & zp = [ 0.0_real64, 1.0_real64, 1.0_real64, 1.0_real64, 2.0_real64, 2.0_real64 ]

! ======================================================================

    iat = nint(a)
    enm = enmax(iat-1)
    scf = scalf(iat-1)
    if (j > 6) then
       write(gsmObj%io%message, 1100)"42"
       call gsmObj%io%print(3, 3, gsmObj%io%message)
    end if
    if (j == 1 .and. e < enm .and. (iat > 2 .and. iat < 300)) then
!   Kalbach's renormalized parametrization is used only
!   for neutron at e<enm
       sinvla = gsmObj%cskalb (1, zro, one, z, a, e)*scf
    else
!   NASA's parametrization is used for p,d,t,he-3,he-4 and
!   for neutron at e>enm
       ea = e/ap(j)
       sinvla = gsmObj%xabs (zp(j), ap(j), z, a, ea)
    endif
    return

! ======================================================================
1100 format("Exceeded ap() array in 'siginv.f90', line ", A)
! ======================================================================
  end function sinvla

  function xabs (gsmObj, zp, ap, zt, at, e)

! ======================================================================
!
!   Absorption xsec revised version RKT-97/4;
!   neutron data from Barashenkov.
!   This gives absorption xsec for given zp, ap, zt, at, e
!   (Mev/nucleon); it can be used for neutrons also.
!
!> Link to the citation:
!> https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19980003828.pdf
!> \todo Made a *.bib file and add a \c \cite command to here.
!
!     Modified by A. J. Sierk, LANL T-16, June, 2002.
!     Modified by A. J. Sierk, October, 2003.
!     Edited by AJS, LANL T-2, December, 2011.
!     Edited by LMK, XCP-3, July 2013 (included error protection)
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use gsm_params, only: zro, one, two, thr, four, pi, twpi, ato3rd

    implicit none
    class(GSM),      intent(inout) :: gsmObj
    real(real64),    intent(in   ) :: zp
    real(real64),    intent(in   ) :: ap
    real(real64),    intent(in   ) :: zt
    real(real64),    intent(in   ) :: at
    real(real64),    intent(in   ) :: e
    real(real64)                   :: xabs

    integer(int32) :: iap, iat, izp, izt
    real(real64)   :: abi, apthrd, atthrd, aux1, aux2, axe, bcm, beta, &
         & bigb, bigr, ce, const, const1, delta, density, dumm1, e20,  &
         & e40, ecm, ecm13, ecmp, ecmt, gcm, plab, rp, rt, sl, t1, temp, &
         & term1, twxsec, vp, vt, x1, xat, xm, xq1, xq2, xq3, xq4, xzt

! ======================================================================

    real(real64), parameter :: c13 = one/thr
    real(real64), parameter :: c43 = four*pi/thr

! ======================================================================

!   Save original particle identities:
    xq1 = zp
    xq2 = ap
    xq3 = zt
    xq4 = at
!   Reverse identities so smaller Z particle is the projectile:
    if (xq3 < xq1) then
       xzt = xq1
       xat = xq2
       xq1 = xq3
       xq2 = xq4
       xq3 = xzt
       xq4 = xat
    endif
    iap = nint(xq2)
    iat = nint(xq4)
    izp = nint(xq1)
    izt = nint(xq3)
!
!  Nucleon-nucleon inelastic xsec not included here:
!
    if (izp*izt == 1 .or. (izp + izt) == 1) then
       xabs = zro
       return
    endif
    rp = gsmObj%radius (iap)
    rt = gsmObj%radius (iat)
    vp =  c43*rp**3
    vt =  c43*rt**3
    if (vp < div0Lim .and. vp > -div0Lim) then
       vp = div0Lim
       write(gsmObj%io%message,1000) "132"
       call gsmObj%io%print(4, 3, gsmObj%io%message)
    end if
    if (vt < div0Lim .and. vt > -div0Lim) then
       vt = div0Lim
       write(gsmObj%io%message,1000) "133"
       call gsmObj%io%print(4, 3, gsmObj%io%message)
    end if
    density = 0.5d0*((xq2/vp) + (xq4/vt))
    const = 1.75d0*density/8.824728d-2
    if (xq4 < div0Lim .and. xq4 > -div0Lim) then
       xq4 = div0Lim
       write(gsmObj%io%message,1000) "139"
       call gsmObj%io%print(4, 3, gsmObj%io%message)
    end if
    abi = one - two*xq3/xq4
    apthrd = ato3rd(iap)
    atthrd = ato3rd(iat)
    if (iap == 1) const = 2.05d0
    if (izp == 2 .and. iap == 4) then
       const1 = 2.77d0 - xq4*8.0d-3 + xq4*xq4*1.8d-5
    endif
    if (izp == 3) const = const/thr
    t1 = 40.d0
    axe = one + e/938.d0
    temp = xq2**2 + xq4**2 + two*xq2*axe*xq4
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(gsmObj%io%message,1000) "154"
       call gsmObj%io%print(4, 3, gsmObj%io%message)
    end if
    gcm = (xq2*axe + xq4)/sqrt(abs(temp))
    if (gcm < div0Lim .and. gcm > -div0Lim) then
       gcm = div0Lim
       write(gsmObj%io%message,1000) "159"
       call gsmObj%io%print(4, 3, gsmObj%io%message)
    end if
    bcm = sqrt(abs(one - one/gcm**2))
    plab = xq2*sqrt(abs(two*938.d0*e + e*e))
    ecmp = gcm*axe*938.d0*xq2 - bcm*gcm*plab - xq2*938.d0
    ecmt = 938.d0*xq4*(gcm - one)
    ecm = ecmp + ecmt
    ecm13 = ecm**c13
    if (ecm13 < div0Lim .and. ecm13 > -div0Lim) then
       ecm13 = div0Lim
       write(gsmObj%io%message,1000) "169"
       call gsmObj%io%print(4, 3, gsmObj%io%message)
    end if
    bigr = rp + rt + 1.2d0*(apthrd + atthrd)/ecm13
    if (bigr < div0Lim .and. bigr > -div0Lim) then
       bigr = div0Lim
       write(gsmObj%io%message,1000) "174"
       call gsmObj%io%print(4, 3, gsmObj%io%message)
    end if
    bigb = 1.43996517d0*xq1*xq3/bigr
    if (izp == 1 .and. iat > 56) bigb = 0.90d0*bigb
    if (iap == 1) then
!   Proton:  (neutron has bigb = 0.)
       if (iat < 4) then
          bigb = 21.d0*bigb
       elseif (iat == 4) then
          bigb = 27.d0*bigb
       elseif (iat == 12) then
          bigb = 3.5d0*bigb
       elseif (iat <= 16 .and. iat >= 13) then
          bigb = (xq4/7.d0)*bigb
       endif
       if (izt == 12) then
          bigb = 1.8d0*bigb
       elseif (izt == 14) then
          bigb = 1.4d0*bigb
       elseif (izt == 20) then
          bigb = 1.3d0*bigb
       endif
    endif
!  Possible for d or t incident on H
    if (iap < 4 .and. iat == 1) bigb = 21.d0*bigb
    xm = one
    if (izp == 0) then
!  Neutron:
       if (iat >= 11 .and. iat < 40) t1 = 30.d0
       if (izt == 14) t1 = 35.d0
       if (izt == 26) t1 = 30.d0
       if (izt == 0) then
          bigb = zro
       elseif (iat < 200) then
          x1 = 2.83d0 - 3.1d-2*xq4 + 1.7d-4*xq4*xq4
          x1 = max (one, x1)
          if (iat < 12) then
             sl = 0.6d0
          elseif (iat == 12) then
             sl = 1.6d0
          else
             sl = one
          endif
          temp = sl*x1
          if (temp < div0Lim .and. temp > -div0Lim) then
             temp = div0Lim
             write(gsmObj%io%message,1000) "220"
             call gsmObj%io%print(4, 3, gsmObj%io%message)
          end if
          xm = (one - x1*exp(-e/(temp)))
       else
!   7/11/96
          xm = (  one - 0.3d0 * exp( -(e - one)/15.d0 )  ) * &
               & (  one - exp( -(e - 0.9d0) )  )
       endif
    endif
    temp = one + exp((250.d0 - e)/75.d0)
    aux1 = one/(temp)
    if (izp == 2 .and. iap == 4) then
!  Alpha:
       const = const1 - 0.8d0*aux1
    endif
    if (izp == 1 .and. iap == 1) then
!  Proton:
       const = 2.05d0 - 0.05d0*aux1
       if (iat < 4) then
          t1 = 55.d0
          const = 1.7d0
       elseif (iat > 45) then
          t1 = 40.d0 + xq4/thr
       endif
       e20 = exp((20.d0 - e)/10.d0)
       temp = one + e20
       aux2 = e20/(temp)
       if (izt == 12) then
          t1 = 40.d0
          const = 2.05d0 - thr*aux2
       elseif (izt == 14) then
          t1 = 40.d0
          const = 2.05d0 - 1.75d0*aux2
       elseif (izt == 18) then
          t1 = 40.d0
          const = 2.05d0 - two*aux2
       elseif (izt == 20) then
          t1 = 40.d0
          e40 = exp((40.d0 - e)/10.d0)
          temp = one + e40
          const = 2.05d0 - one*e40/(temp)
       endif
    elseif (izp == 0 .and. iap == 1) then
!  Neutron:
       if (density < div0Lim .and. density > -div0Lim) then
          density = div0Lim
          write(gsmObj%io%message,1000) "278"
          call gsmObj%io%print(4, 3, gsmObj%io%message)
       end if
       const = two*0.134457d0/density
       if (iat < 60) then
          const = const - 1.5d0*abi
          ! Adjust cross section for shell corrections (NOTE for large E it's
          ! negligible)
          if (iat <= 40 .and. ((0.01_real64 * e - 1.7_real64) <= 12)) then
             temp = one + exp(-1.70d0 + 0.01d0*e)
             const = const + 0.25d0/(temp)
         end if
       elseif (iat > 140 .and. iat < 200) then
          const = const - 1.5d0*abi
       endif
       temp = xq4 - xq3
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(gsmObj%io%message,1000) "295"
          call gsmObj%io%print(4, 3, gsmObj%io%message)
       end if
       if (izt > 82) const = const - xq3/(temp)
       if (izt <= 20 .and. izt >= 10) then
          if ((0.1_real64 * e - two) <= 12) then
             temp = one + exp(0.1d0*e - two)
             const = const - one/(temp)
          end if
       elseif (izt >= 82) then
          temp = one + exp(0.05d0*e - one)
          const = const - two/(temp)
       endif
    endif
    if (t1 < div0Lim .and. t1 > -div0Lim) then
       t1 = div0Lim
       write(gsmObj%io%message,1000) "316"
       call gsmObj%io%print(4, 3, gsmObj%io%message)
    end if
    ce = const*(one - exp(-e/t1)) - 0.292d0*exp(-e/792.d0) * &
         & cos(0.229d0*e**0.453d0)
    dumm1 = atthrd + apthrd
    if (dumm1 < div0Lim .and. dumm1 > -div0Lim) then
       dumm1 = div0Lim
       write(gsmObj%io%message,1000) "323"
       call gsmObj%io%print(4, 3, gsmObj%io%message)
    end if
    term1 = atthrd*apthrd/dumm1
    delta = 1.615d0*term1 - 0.873d0*ce
    if (ecm13 < div0Lim .and. ecm13 > -div0Lim) then
       ecm13 = div0Lim
       write(gsmObj%io%message,1000) "329"
       call gsmObj%io%print(4, 3, gsmObj%io%message)
    end if
    delta = delta + 0.140d0*term1/ecm13
    if (xq2 < div0Lim .and. xq2 > -div0Lim) then
       xq2 = div0Lim
       write(gsmObj%io%message,1000) "334"
       call gsmObj%io%print(4, 3, gsmObj%io%message)
    end if
    delta = delta + 0.794d0*abi*xq1/xq2
    delta = -delta
    beta = one
    twxsec = 10.d0*pi*1.26d0*1.26d0*beta*(0.873d0*dumm1 - delta)**2
    if (ecm < div0Lim .and. ecm > -div0Lim) then
       ecm = div0Lim
       write(gsmObj%io%message,1000) "342"
       call gsmObj%io%print(4, 3, gsmObj%io%message)
    end if
    xabs = twxsec*(one - bigb/ecm)*xm
    xabs = max (xabs, zro)

    return

! ======================================================================
1000 format("Divide by zero error prevented in 'siginv.f90', line(s) ", A)
! ======================================================================
  end function xabs

  function radius (gsmObj, ia)

! ======================================================================
!
!     Edited in September, 2003 by A. J. Sierk
!     Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use gsm_params, only: ato3rd

    implicit none
    class(GSM),     intent(inout) :: gsmObj
    integer(int32), intent(in   ) :: ia
    real(real64)                  :: radius

    integer(int32) :: i

! ======================================================================

    integer(int32), parameter, dimension(23) :: na =             &
         & [ 1,     2,     3,     4,     6,     7,     9,    10, &
         &  11,    12,    13,    14,    15,    16,    17,    18, &
         &  19,    20,    22,    23,    24,    25,    26         ]

    real(real64),   parameter, dimension(23) :: rms =              &
         & [ 0.850d0, 2.095d0, 1.976d0, 1.671d0, 2.570d0, 2.410d0, &
         &   2.519d0, 2.45d0 , 2.420d0, 2.471d0, 2.440d0, 2.580d0, &
         &   2.611d0, 2.730d0, 2.662d0, 2.727d0, 2.900d0, 3.040d0, &
         &   2.969d0, 2.940d0, 3.075d0, 3.110d0, 3.060d0           ]

    real(real64), parameter :: fact = sqrt( 5.0_real64 / 3.0_real64 )

! ======================================================================

    if (ia == 5 .or. ia == 8 .or. ia == 21 .or. ia > 26) then
       radius = fact * (0.84d0*ato3rd(ia)  + 0.55d0)
    else
       do i = 1,23
          if (ia == na(i)) radius = fact*rms(i)
       end do
    endif

    return

! ======================================================================
  end function radius


  function cskalb (gsmObj, kp, zp, ap, zt, at, elab)

! ======================================================================
!
!     Written in 1982; revised 1990
!
!     Calculate optical model reaction cross sections
!     using empirical parameterization
!     of Narasimha Murthy, Chaterjee, and Gupta
!     going over to the geometrical limit at high energy
!
!           Proton cross sections scaled down with signor for a<100
!           (appropriate for Becchetti-Greenlees potential)
!           Neutron cross sections scaled down sith signor for a<40
!           (appropriate for Mani et al potential)
!
!     parameter values set in subroutine sigpar
!
!     Edited in September, 2003 by A. J. Sierk
!     Edited by AJS, LANL T-2, December, 2011.
!     Edited by LMK, XCP-3, July 2013 (included error protection)
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use gsm_params, only: zro, one, two, ato3rd
    use generalizedSpallationData, only: xl0, xl1, xm0, xm1, xn0, xn1, &
         & xn2, xp0, xp1, xp2

    implicit none
    class(GSM),     intent(inout) :: gsmObj
    integer(int32), intent(in   ) :: kp
    real(real64),   intent(in   ) :: zp
    real(real64),   intent(in   ) :: ap
    real(real64),   intent(in   ) :: zt
    real(real64),   intent(in   ) :: at
    real(real64),   intent(in   ) :: elab
    real(real64)                  :: cskalb

    integer(int32) :: iat, jout
    real(real64)   :: a, ares, athrd, b, c, cut, ec, ecsq, ecut, ecut2, &
         & etest, flow, geom, p, ra, sig, signor, signor2, spill, &
         & temp, w, xlamb, xmu, xnu, xnulam

! ======================================================================

    flow = 1.d-250
    spill = 1.d+250
    jout = nint (ap)
    iat = nint (at)
    ares = at + ap
    athrd = ato3rd(iat + jout)
    signor = one
!    signor reduces p and n result for light targs as per expt.
    if (kp == 1) then
       if (ares < 40.d0) signor = 0.7d0 + ares*0.0075d0
       xlamb = xl0(1)/athrd + xl1(1)
       xmu = xm0(1)*athrd + xm1(1)*athrd*athrd
       xnu = xn0(1)*athrd*ares + xn1(1)*athrd*athrd + xn2(1)
       ec = 0.5d0
       ecsq = 0.25d0
       p = xp0(1)
       xnulam = one
       etest = 32.d0
!     etest is the energy above which the rxn cross section is
!     compared with the geometrical limit and the max taken.
!     xnulam here is a dummy value to be used later.
       ra = zro
    else
       ra = 1.20d0
       if (kp == 2) then
          ra = zro
          if (ares < 60.d0) then
             signor = 0.92d0
          elseif (ares < 100.d0) then
             signor = 0.8d0 + ares*0.002d0
          endif
       endif
       temp = 1.5d0*athrd + ra
       ec = 1.439965173d0*zp*zt/(temp)
       if (ec < div0Lim .and. ec > -div0Lim) then
          ec = div0Lim
          write(gsmObj%io%message,1000) "493"
          call gsmObj%io%print(4, 3, gsmObj%io%message)
       end if
       ecsq = ec * ec
       p = xp0(kp) + xp1(kp)/ec + xp2(kp)/ecsq
       xlamb = xl0(kp)*ares + xl1(kp)
       a = ares**xm1(kp)
       xmu = xm0(kp) * a
       xnu = a*(xn0(kp) + xn1(kp)*ec + xn2(kp)*ecsq)
       if (jout == 2 .or. jout == 3) ra = 0.8d0
!    new values of ra are for calculating the geometrical limit
!    to the cross section.
       if (kp == 2) then
          c = min(3.15d0, ec*0.5d0)
          w = 0.7d0*c/3.15d0
!    c and w are for the global corr'n factor for elab<ec
!    for light targs they are scaled down from global values
       endif
       if (xlamb < div0Lim .and. xlamb > -div0Lim) then
          xlamb = div0Lim
          write(gsmObj%io%message,1000) "511"
          call gsmObj%io%print(4, 3, gsmObj%io%message)
       end if
       xnulam = xnu/xlamb
       if (xnulam > spill) xnulam = zro
       if (xnulam >= flow) then
          if (kp == 2) then
             etest = sqrt(xnulam) + 7.d0
          else
             etest = 1.2d0*sqrt(xnulam)
          endif
! for xnulam > 0, sig reaches a maximum at sqrt(xnulam).
       endif
    endif
    a = -two*p*ec + xlamb - xnu/ecsq
    b = p*ecsq + xmu + two*xnu/ec
    ecut = zro
    cut = a*a - 4.d0*p*b
    if (cut > zro) ecut = sqrt(abs(cut))
    ecut = (ecut - a)/(p + p)
    ecut2 = ecut
    if (cut < zro) ecut2 = ecut - two
    sig = zro
    if (elab <= ec) then
       if (elab > ecut2) then
          sig = (p*elab*elab + a*elab + b) * signor
          if (kp == 2) then
             if (w < div0Lim .and. w > -div0Lim) then
                w = div0Lim
                write(gsmObj%io%message,1000) "539"
                call gsmObj%io%print(4, 3, gsmObj%io%message)
             end if
             signor2 = (ec - elab - c)/w
             signor2 = one + exp(signor2)
             sig = sig/signor2
          endif
!    first signor gives empirical global corr'ns at low elab
!    second signor corrects values near elab=0; light nuclei
       endif
    else
       sig = (xlamb*elab + xmu + xnu/elab) * signor
       geom = zro
       if (xnulam >= flow .and. elab >= etest) then
          geom = sqrt(abs(ap*elab))
          if (geom < div0Lim .and. geom > -div0Lim) then
             geom = div0Lim
             write(gsmObj%io%message,1000) "563"
             call gsmObj%io%print(4, 3, gsmObj%io%message)
          end if
          geom = 1.23d0*athrd + ra + 4.573d0/geom
          geom = 31.416d0*geom*geom
          sig = max(geom, sig)
       endif
    endif

    cskalb = sig
    return
! ======================================================================
1000 format("Divide by zero error prevented in 'siginv.f90', line(s) ", A)
! ======================================================================
  end function cskalb
