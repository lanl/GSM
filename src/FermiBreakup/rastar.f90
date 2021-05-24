
  subroutine rastar (fbuObj, ap, zp, up, pn, nst, nw, results)

! ======================================================================
!
!    Called by: FERMIDEC
!
!    Calls: CLPV, DISIMP, PINT, RAZVAL, WEPT
!
!    Last change: 13-Aug-2003 BY NVMokhov
!    Edited by A. J. Sierk, LANL T-16, September, 2003.
!    Modified by K. K. Gudima, November, 2004.
!    Edited by AJS, January, 2005.
!    Edited by A. J. Sierk, LANL T-2, February, 2009.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection)
!    Edited by CMJ, XCP-3, July 2018 (creation of FermiBreakup class).
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use fermiBreakupParams, only : zro, thousand, thousandth, &
         & electron_mass, nucleon_mass, proton_mass
    use fermiBreakUpData, only : dm, min_fermi_AData

    implicit none
    class(FermiBreakup), intent(inout) :: fbuObj
    real(real64),        intent(in   ) :: ap
    real(real64),        intent(in   ) :: zp
    real(real64),        intent(in   ) :: up
    real(real64),        intent(in   ), dimension(3) :: pn
    integer(int32),      intent(inout) :: nst ! # Stable fragments produced from FBU process
    integer(int32),      intent(in   ) :: nw
    type(fermiBreakUpResults), intent(inout) :: results

    integer(int32) :: i, ia, inp, iz, izp, j, jaz, k, l, &
         & lv, m1, msa, mv, n, nb, nf, ng, npt, nv
    real(real64)   :: br, cf, cmn, cmp, ct, daz, eap, en, pq, sf, &
         & sr, st, tk, tmn, tn

    real(real64), dimension(3) :: vn, pnc, pnl
    real(real64), dimension(3, fbuObj%options%recNumNucleons) :: pn0
    integer(int32), dimension(fbuObj%options%recNumNucleons) :: ipa, ipz, nc
    real(real64), dimension(fbuObj%options%recNumNucleons) :: amk

    type(fbuInterimCalc) :: calc

! ======================================================================

    izp = nint(zp) + 1
    inp = nint(ap - zp) + 1
    !   KKG 11/02/04:
    if (inp <= 0 .or. izp <= 0 .or. inp > min_fermi_AData .or. izp > (min_fermi_AData-1) ) then
       ! Particle index of the 'dm(X,Y)' array bounds
       cmp = proton_mass*ap
    else
       cmp = nucleon_mass*ap + thousandth*dm(inp,izp) - electron_mass*zp
    endif

    ! Setup momenta, nucleus veloctiy
    pq = pn(1)**2 + pn(2)**2 + pn(3)**2
    eap = sqrt(pq + cmp**2)
    if (eap < div0Lim .and. eap > -div0Lim) then
       eap = div0Lim
       write(fbuObj%io%message,1000) "56"
       call fbuObj%io%print(4, 3, fbuObj%io%message)
    end if
    do k = 1,3
       vn(k) = pn(k)/eap
    end do


    ! Start of calculation
    if (ap >= 1.9d0) then
       call fbuObj%razval (ap, zp, up, nf, npt, calc)
       if (nf >= 1) then
          ia = nint(ap)
          iz = nint(zp)
          jaz = 100*ia + iz
          daz = fbuObj%wept (ia, iz) + up
          do l = 1,nw
             br = fbuObj%rang()
             sr = zro
             do n = 1,nf
                nb = n
                sr = sr + calc%wnq(n)
                if (br < sr) go to 10
             end do
10           nc(l) = nb
          end do
          nv = nw
          k = 0
          msa = 0
          ng = 0
          do n = 1,npt
             k = k + 1
             ipa(k) = calc%iaz(n)/100
             ipz(k) = calc%iaz(n) - 100*ipa(k)
             msa = msa + calc%iaz(n)
             if (msa >= jaz) then
                ng = ng + 1
20              if (nv >= 1) then
                   do l = 1,nv
                      lv = l
                      if (nc(l) == ng) go to 30
                   end do
                   go to 40
30                 nv = nv - 1
                   if (nv >= 1) then
                      do  l = lv,nv
                         nc(l) = nc(l+1)
                      end do
                   endif
                   tn = daz
                   do i = 1,k
                      mv = nst + i
                      if (mv > fbuObj%options%recNumNucleons) then
                         ! Particle will be tracked, but cannot be stored
                         write(fbuObj%io%message, 100) fbuObj%options%recNumNucleons, (nst+k)
                         call fbuObj%io%print(3, 3, fbuObj%io%message)
                         results%simState = results%simState + exceededProgenyArray
                      end if
                      amk(i) = fbuObj%wept (ipa(i), ipz(i))
                      tn = tn - amk(i)
                   end do
                   call fbuObj%disimp (k, amk, pn0, tn)
                   do i = 1,k
                      do j = 1,3
                         pnc(j) = pn0(j,i)
                      end do
                      en = sqrt(pnc(1)**2 + pnc(2)**2 + pnc(3)**2 + amk(i)**2)
                      call fbuObj%clpv (pnc, vn, pnl, en)
                      call fbuObj%pint (pnl, ct, st, cf, sf, tk, amk(i))
                      mv = nst + i

                      ! If stored all possible progeny, return (no more can be stored)
                      if ( mv > results%maxProgeny ) then
                         write(fbuObj%io%message,2000)
                         call fbuObj%io%print(3, 3, fbuObj%io%message)
                         write(fbuObj%io%message,2010) (mv - results%maxProgeny)
                         call fbuObj%io%print(3, 3, fbuObj%io%message)
                         results%simState = results%simState + exceededProgenyArray
                         nst = mv - 1
                         return
                      end if
                      results%progenyBnk(mv)%numBaryons = dble(ipa(i))
                      results%progenyBnk(mv)%numProtons = dble(ipz(i))
                      results%progenyBnk(mv)%linearXMom = pnl(1)*thousand
                      results%progenyBnk(mv)%linearYMom = pnl(2)*thousand
                      results%progenyBnk(mv)%linearZMom = pnl(3)*thousand
                      results%progenyBnk(mv)%kinEnergy  = tk*thousand
                      results%progenyBnk(mv)%restMass   = amk(i)*thousand
                   end do
                   nst = mv
                   go to 20
                endif
40              k = 0
                msa = 0
             endif
          end do
          return
       endif
    endif

    ! Particle undergoing fermi breakup is a nucleon (only 1 progeny, namely itself)
    m1 = nst + 1
    nst = m1
    if (ap > 0.1d0) then
       if ( inp > min_fermi_AData .or. izp > (min_fermi_AData-1) ) then
          ! Array dm(X,Y) will be exceeded; warn client and approximate
          write(fbuObj%io%message,1100) inp, izp
          call fbuObj%io%print(3, 3, fbuObj%io%message)
          write(fbuObj%io%message,1110)
          call fbuObj%io%print(3, 3, fbuObj%io%message)
          if ( inp > min_fermi_AData ) inp = min_fermi_AData           ! Approximate neutron number
          if ( izp > min_fermi_AData - 1 ) izp = min_fermi_AData - 1   ! Approximate proton  number
          results%simState = results%simState + dataApproximated
       end if
       cmn = nucleon_mass*ap + thousandth*dm(inp,izp) - electron_mass*zp
       call fbuObj%pint (pn, ct, st, cf, sf, tmn, cmn)
       results%progenyBnk(m1)%numBaryons = ap
       results%progenyBnk(m1)%numProtons = zp
       results%progenyBnk(m1)%linearXMom = pn(1)*thousand
       results%progenyBnk(m1)%linearYMom = pn(2)*thousand
       results%progenyBnk(m1)%linearZMom = pn(3)*thousand
       results%progenyBnk(m1)%kinEnergy  = tmn*thousand
       results%progenyBnk(m1)%restMass   = cmn*thousand      
    endif

    return
! ======================================================================
100 format ("Only ", i5, " of the ", i5, " Fermi Break-Up progeny ", &
         & "can be tracked (not enough memory).")
1000 format("Divide by zero error in 'rastar.f90', line ", A)
1100 format("The Fermi Break-Up model's 'dm' data array was ", &
          & "exceeded (N=", i3, ", Z=", i3, ").")
1110 format("   Using closest valid value.")
2000 format("The internal fermi break-up fragment array size is exceeded.")
2010 format("   The remaining ", i3, " fragment(s) will not be tallied.")
! ======================================================================
  end subroutine rastar
