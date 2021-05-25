
  subroutine prtdist (gsmObj, output, ipar, om, const, doo, siginn)

! ======================================================================
!
!   This subroutine prints out spectra, angular distriputions, and
!   double differential cross sections for particle of type ipar.
!
!   Extracted from old TYPEOUT by A. J. Sierk, LANL T-16, October, 2003.
!   Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use gsm_params, only: zro, one

    implicit none
    class(GSM),     intent(inout) :: gsmObj
    class(GSMOutput), intent(in   ) :: output
    integer(int32), intent(in   ) :: ipar
    real(real64),   intent(in   ) :: om(20)
    real(real64),   intent(in   ) :: const(20)
    real(real64),   intent(in   ) :: doo(10)
    real(real64),   intent(in   ) :: siginn

    integer(int32) :: i, i1, i2, i3, i4, i5, j, k
    real(real64)   :: a1, a2, a3, a4, a5a, a6a, a7a, an1, da1, da2, &
         & da3, da4, da5a, da6a, da7a, ds1, ds1n, ds2, ds2n, ds3, &
         & ds3n, ds4, ds4n, ds5a, ds5n, ds6a, ds6n, ds7a, ds7n, s1, &
         & s1n, s2, s2n, s3, s3n, s4, s4n, s5a, s5n, s6a, s6n, s7a, s7n
    logical :: ls1, ls2, prnt1d, prnt2d

    real(real64), dimension(200) :: s5=zro, ds5=zro, s6=zro, ds6=zro, &
         & s7=zro, ds7=zro, a5=zro, da5=zro, a6=zro, da6=zro, a7=zro, &
         & da7=zro

! ======================================================================

    character(LEN=3), parameter, dimension(9) :: par = &
         & ['n  ', 'p  ', 'd  ', 't  ', 'He3', 'He4', 'pi-', 'pi0', &
         &  'pi+' ]

    integer(int32), parameter, dimension(9, 4) :: ipp = reshape( &
         & ([  1,  2,  3,  1, &
         &     4,  5,  6,  2, &
         &     7,  8,  8,  3, &
         &     9, 10, 10,  4, &
         &    11, 12, 12,  5, &
         &    13, 14, 14,  6, &
         &     7,  7,  7,  7, &
         &     8,  8,  8,  8, &
         &     9,  9,  9,  9  ] ), shape(ipp), order=([2, 1])  )

    real(real64)   :: te1, te2, dtt, se, dtdo, d2spec, d2spe
    integer(int32) ::  nt2, nt3, ntet, nti, ntt
    common /d2sdto/  te1(200), te2(200), dtt(200), &
         & se(200), dtdo(200,10), d2spec(9,10,200), &
         & d2spe(14,10,200), ntt, ntet, nt2, nt3, nti(4)
    real(real64)   :: speco, angco, d2speco
    common /rescoa/  speco(4,200), angco(4,20), d2speco(4,10,200)
    real(real64)   :: spef, angf, d2spef
    common /resfis/  spef(6,200), angf(6,20), d2spef(6,10,200)
    real(real64)   :: spe, an
    common /resmac/  spe(14,200), an(14,20)
    real(real64)   :: specsp, angsp, d2spesp
    common /respal/  specsp(6,200), angsp(6,20), d2spesp(6,10,200)
    real(real64)   :: specpf, angpf, d2spepf
    common /resprf/  specpf(6,200), angpf(6,20), d2spepf(6,10,200)
    real(real64)   :: spec, ang, chan, dadz
    common /result/  spec(9,200), ang(9,20), chan(218), dadz(351,151)

! ======================================================================

    i1 = ipp(ipar,1)
    i2 = ipp(ipar,2)
    i3 = ipp(ipar,3)
    i4 = ipp(ipar,4)
    i5 = i4 - 2
    if (output%energySpectra > 0) then
       write (31, 1900)
       prnt1d = (dble(int(dtt(1)*10 + 1.d-6))/10.d0)  ==  dtt(1)
       prnt2d = .not.prnt1d
       if (prnt1d) write (31, 2100) te1(1), te2(ntt)
       if (prnt2d) write (31, 2200) te1(1), te2(ntt)
       if (ipar < 3) then
          write (31, 1000) par(ipar)
       elseif (ipar >= 3 .and. ipar < 7) then
          write (31, 1700) par(ipar)
       elseif (ipar > 6) then
          write (31, 1400) par(ipar)
       endif
       do k = 1,ntt
          prnt1d = (dble(int(dtt(k)*10 + 1.d-6))/10.d0)  ==  dtt(k)
          prnt2d = .not.prnt1d
          s1 = spec(i4,k)*se(k)
          ds1 = sqrt(abs(spec(i4,k)))*se(k)
          if (ipar <= 2) then
             s2 = spe(i1,k)*se(k)
             ds2 = sqrt(abs(spe(i1,k)))*se(k)
             s3 = spe(i2,k)*se(k)
             ds3 = sqrt(abs(spe(i2,k)))*se(k)
             s4 = spe(i3,k)*se(k)
             ds4 = sqrt(abs(spe(i3,k)))*se(k)
          elseif (ipar > 2 .and. ipar <= 6) then
             s2 = speco(i5,k)*se(k)
             ds2 = sqrt(abs(speco(i5,k)))*se(k)
             s3 = spe(i1,k)*se(k)
             ds3 = sqrt(abs(spe(i1,k)))*se(k)
             s4 = spe(i2,k)*se(k)
             ds4 = sqrt(abs(spe(i2,k)))*se(k)
          endif
          if (ipar < 7) then
             s5(k) = specsp(i4,k)*se(k)
             ds5(k) = sqrt(abs(specsp(i4,k)))*se(k)
             s6(k) = specpf(i4,k)*se(k)
             ds6(k) = sqrt(abs(specpf(i4,k)))*se(k)
             s7(k) = spef(i4,k)*se(k)
             ds7(k) = sqrt(abs(spef(i4,k)))*se(k)
             ls1 = s1 > zro .or. s2 > zro .or. s3 > zro .or. &
                  & s4 > zro
             ls2 = s5(k) > zro .or. s6(k) > zro .or. s7(k) > zro
             if (ls1) then
                if (prnt1d) write (31, 2400) te1(k), te2(k), s1, &
                     & ds1, s2, ds2, s3, ds3, s4, &
                     & ds4
                if (prnt2d) write (31, 2500) te1(k), te2(k), s1, &
                     & ds1, s2, ds2, s3, ds3, s4, &
                     & ds4
             endif
          else
             if (s1 > zro) then
                if (prnt1d) write (31, 2400) te1(k), te2(k), s1, ds1
                if (prnt2d) write (31, 2500) te1(k), te2(k), s1, ds1
             endif
          endif
       end do
       s1 = spec(i4,nt3)*siginn
       ds1 = sqrt(abs(spec(i4,nt3)))*siginn
       if (ipar < 7) then
          if (ipar >= 3) then
             s2 = speco(i5,nt3)*siginn
             ds2 = sqrt(abs(speco(i5,nt3)))*siginn
             s3 = spe(i1,nt3)*siginn
             ds3 = sqrt(abs(spe(i1,nt3)))*siginn
             s4 = spe(i2,nt3)*siginn
             ds4 = sqrt(abs(spe(i2,nt3)))*siginn
          elseif (ipar < 3) then
             s2 = spe(i1,nt3)*siginn
             ds2 = sqrt(abs(spe(i1,nt3)))*siginn
             s3 = spe(i2,nt3)*siginn
             ds3 = sqrt(abs(spe(i2,nt3)))*siginn
             s4 = spe(i3,nt3)*siginn
             ds4 = sqrt(abs(spe(i3,nt3)))*siginn
          endif
          write (31, 3000) s1, ds1, s2, ds2, s3, ds3, s4, ds4
       elseif (ipar >= 7) then
          write (31, 3000) s1, ds1
       endif
       if (output%energySpectra == 2 .and. ipar < 7) then
          write (31, 1100) par(ipar)
          do k = 1,ntt
             prnt1d = (dble(int(dtt(k)*10 + 1.d-6))/10.d0)  ==  dtt(k)
             prnt2d = .not.prnt1d
             ls2 = s5(k) > zro .or. s6(k) > zro .or. s7(k) > zro
             if (ls2) then
                if (prnt1d) &
                     & write (31, 2400) te1(k), te2(k), s5(k), ds5(k), s6(k), &
                     & ds6(k), s7(k), ds7(k)
                if (prnt2d) &
                     & write (31, 2500) te1(k), te2(k), s5(k), ds5(k), s6(k), &
                     & ds6(k), s7(k), ds7(k)
             endif
          end do
          s5a = specsp(i4,nt3)*siginn
          ds5a = sqrt(abs(specsp(i4,nt3)))*siginn
          s6a = specpf(i4,nt3)*siginn
          ds6a = sqrt(abs(specpf(i4,nt3)))*siginn
          s7a = spef(i4,nt3)*siginn
          ds7a = sqrt(abs(spef(i4,nt3)))*siginn
          write (31, 3000) s5a, ds5a, s6a, ds6a, s7a, ds7a
       endif
       if (ipar == 1) then
!   Print out normalized neutron spectrum components: (AJS, 9/26/03)
          if ( s1 < div0Lim .and. s1 > -div0Lim ) then
             s1 = div0Lim
             write(gsmObj%io%message, 5000) "188"
             call gsmObj%io%print(4, 3, gsmObj%io%message)
          end if
          an1 = one/s1
          write (31, 2000)
          prnt1d = (dble(int(dtt(1)*10 + 1.d-6))/10.d0)  ==  dtt(1)
          prnt2d = .not.prnt1d
          if (prnt1d) write (31, 2100) te1(1), te2(ntt)
          if (prnt2d) write (31, 2200) te1(1), te2(ntt)
          write (31, 1000) par(ipar)
          do k = 1,ntt
             prnt1d = (dble(int(dtt(k)*10 + 1.d-6))/10.d0)  ==  dtt(k)
             prnt2d = .not.prnt1d
             s1n = spec(i4,k)*se(k)*an1
             ds1n = sqrt(abs(spec(i4,k)))*se(k)*an1
             s2n = spe(i1,k)*se(k)*an1
             ds2n = sqrt(abs(spe(i1,k)))*se(k)*an1
             s3n = spe(i2,k)*se(k)*an1
             ds3n = sqrt(abs(spe(i2,k)))*se(k)*an1
             s4n = spe(i3,k)*se(k)*an1
             ds4n = sqrt(abs(spe(i3,k)))*se(k)*an1
             ls1 = s1n > zro .or. s2n > zro .or. s3n > zro .or. &
                  & s4n > zro
             if (ls1) then
                if (prnt1d) write (31, 2400) te1(k), te2(k), s1n, ds1n, &
                     & s2n, ds2n, s3n, ds3n, s4n, ds4n
                if (prnt2d) write (31, 2500) te1(k), te2(k), s1n, ds1n, &
                     & s2n, ds2n, s3n, ds3n, s4n, ds4n
             endif
          end do
          s1n = spec(i4,nt3)*siginn*an1
          ds1n = sqrt(abs(spec(i4,nt3)))*siginn*an1
          s2n = spe(i1,nt3)*siginn*an1
          ds2n = sqrt(abs(spe(i1,nt3)))*siginn*an1
          s3n = spe(i2,nt3)*siginn*an1
          ds3n = sqrt(abs(spe(i2,nt3)))*siginn*an1
          s4n = spe(i3,nt3)*siginn*an1
          ds4n = sqrt(abs(spe(i3,nt3)))*siginn*an1
          write (31, 3000) s1n, ds1n, s2n, ds2n, s3n, ds3n, s4n, ds4n
          if (output%energySpectra == 2) then
             write (31, 1100) par(ipar)
             do k = 1,ntt
                prnt1d = (dble(int(dtt(k)*10 + 1.d-6))/10.d0)  ==  dtt(k)
                prnt2d = .not.prnt1d
                s5(k) = specsp(i4,k)*se(k)*an1
                ds5(k) = sqrt(abs(specsp(i4,k)))*se(k)*an1
                s6(k) = specpf(i4,k)*se(k)*an1
                ds6(k) = sqrt(abs(specpf(i4,k)))*se(k)*an1
                s7(k) = spef(i4,k)*se(k)*an1
                ds7(k) = sqrt(abs(spef(i4,k)))*se(k)*an1
                ls2 = s5(k) > zro .or. s6(k) > zro .or. s7(k) > zro
                if (ls2) then
                   if (prnt1d) write (31, 2400) te1(k), te2(k), s5(k), &
                        & ds5(k), s6(k), ds6(k), s7(k), ds7(k)
                   if (prnt2d) write (31, 2500) te1(k), te2(k), s5(k), &
                        & ds5(k), s6(k), ds6(k), s7(k), ds7(k)
                endif
             end do
             s5n = specsp(i4,nt3)*siginn*an1
             ds5n = sqrt(abs(specsp(i4,nt3)))*siginn*an1
             s6n = specpf(i4,nt3)*siginn*an1
             ds6n = sqrt(abs(specpf(i4,nt3)))*siginn*an1
             s7n = spef(i4,nt3)*siginn*an1
             ds7n = sqrt(abs(spef(i4,nt3)))*siginn*an1
             write (31, 3000) s5n, ds5n, s6n, ds6n, s7n, ds7n
          endif
!   End of normalized neutron energy spectrum ^
       endif
    endif
!   End of spectrum loop ^
    if (output%angularSpectra > 0) then
       write (31, 2300)
       if (ipar < 3) then
          write (31, 1500) par(ipar)
       elseif (ipar >= 3 .and. ipar < 7) then
          write (31, 1550) par(ipar)
       elseif (ipar >= 7) then
          write (31, 1300) par(ipar)
       endif
       do j = 1,18
          a1 = ang(i4,j)*const(j)
          da1 = sqrt(abs(ang(i4,j)))*const(j)
          if (ipar < 3) then
             a2 = an(i1,j)*const(j)
             da2 = sqrt(abs(an(i1,j)))*const(j)
             a3 = an(i2,j)*const(j)
             da3 = sqrt(abs(an(i2,j)))*const(j)
             a4 = an(i3,j)*const(j)
             da4 = sqrt(abs(an(i3,j)))*const(j)
          elseif (ipar >= 3 .and. ipar < 7) then
             a2 = angco(i5,j)*const(j)
             da2 = sqrt(abs(angco(i5,j)))*const(j)
             a3 = an(i1,j)*const(j)
             da3 = sqrt(abs(an(i1,j)))*const(j)
             a4 = an(i2,j)*const(j)
             da4 = sqrt(abs(an(i2,j)))*const(j)
          elseif (ipar >= 7) then
             if (a1 > zro) write (31, 2700) om(j), a1, da1
          endif
          if (ipar < 7) then
             ls1 = a1 > zro .or. a2 > zro .or. a3 > zro .or. a4 > zro
             if (ls1) write (31, 2700) om(j), a1, da1, a2, da2, a3, da3, &
                  & a4, da4
          endif
       end do
       a1 = ang(i4,20)*siginn
       da1 = sqrt(abs(ang(i4,20)))*siginn
       if (ipar < 3) then
          a2 = an(i1,20)*siginn
          da2 = sqrt(abs(an(i1,20)))*siginn
          a3 = an(i2,20)*siginn
          da3 = sqrt(abs(an(i2,20)))*siginn
          a4 = an(i3,20)*siginn
          da4 = sqrt(abs(an(i3,20)))*siginn
          write (31, 3200) a1, da1, a2, da2, a3, da3, a4, da4
       elseif (ipar >= 3 .and. ipar < 7) then
          a2 = angco(i5,20)*siginn
          da2 = sqrt(abs(angco(i5,20)))*siginn
          a3 = an(i1,20)*siginn
          da3 = sqrt(abs(an(i1,20)))*siginn
          a4 = an(i2,20)*siginn
          da4 = sqrt(abs(an(i2,20)))*siginn
          write (31, 3200) a1, da1, a2, da2, a3, da3, a4, da4
       elseif (ipar >= 7) then
          write (31, 3200) a1, da1
       endif
       if (output%angularSpectra == 2 .and. ipar < 7) then
          write (31, 1600) par(ipar)
          do j = 1,18
             a5(j) = angsp(i4,j)*const(j)
             da5(j) = sqrt(abs(angsp(i4,j)))*const(j)
             a6(j) = angpf(i4,j)*const(j)
             da6(j) = sqrt(abs(angpf(i4,j)))*const(j)
             a7(j) = angf(i4,j)*const(j)
             da7(j) = sqrt(abs(angf(i4,j)))*const(j)
             ls2 = a5(j) > zro .or. a6(j) > zro .or. a7(j) > zro
             if (ls2) write (31, 2700) om(j), a5(j), da5(j), a6(j), &
                  & da6(j), a7(j), da7(j)
          end do
          a5a = angsp(i4,20)*siginn
          da5a = sqrt(abs(angsp(i4,20)))*siginn
          a6a = angpf(i4,20)*siginn
          da6a = sqrt(abs(angpf(i4,20)))*siginn
          a7a = angf(i4,20)*siginn
          da7a = sqrt(abs(angf(i4,20)))*siginn
          write (31, 3200) a5a, da5a, a6a, da6a, a7a, da7a
       endif
    endif
!   End of angular distribution loop ^
    if (output%doubleDiffSpectra > 0) then
       do i = 1,ntet
          if (d2spec(i4,i,nt3) == zro) return
          write (31, 1800) output%angleBins(i)%lowerBound, output%angleBins(i)%upperBound
          if (ipar < 3) then
             write (31, 1000) par(ipar)
          elseif (ipar >= 3 .and. ipar < 7) then
             write (31, 1700) par(ipar)
          elseif (ipar > 6) then
             write (31, 1400) par(ipar)
          endif
          do k = 1,ntt
             prnt1d = (dble(int(dtt(k)*10 + 1.d-6))/10.d0)  ==  dtt(k)
             prnt2d = .not.prnt1d
             s1 = d2spec(i4,i,k)*dtdo(k,i)
             ds1 = sqrt(abs(d2spec(i4,i,k)))*dtdo(k,i)
             if (ipar < 3) then
                s2 = d2spe(i1,i,k)*dtdo(k,i)
                ds2 = sqrt(abs(d2spe(i1,i,k)))*dtdo(k,i)
                s3 = d2spe(i2,i,k)*dtdo(k,i)
                ds3 = sqrt(abs(d2spe(i2,i,k)))*dtdo(k,i)
                s4 = d2spe(i3,i,k)*dtdo(k,i)
                ds4 = sqrt(abs(d2spe(i3,i,k)))*dtdo(k,i)
             elseif (ipar >= 3 .and. ipar < 7) then
                s2 = d2speco(i5,i,k)*dtdo(k,i)
                ds2 = sqrt(abs(d2speco(i5,i,k)))*dtdo(k,i)
                s3 = d2spe(i1,i,k)*dtdo(k,i)
                ds3 = sqrt(abs(d2spe(i1,i,k)))*dtdo(k,i)
                s4 = d2spe(i2,i,k)*dtdo(k,i)
                ds4 = sqrt(abs(d2spe(i2,i,k)))*dtdo(k,i)
             elseif (ipar >= 6) then
                if (s1 > zro) then
                   if (prnt1d) write (31, 2400) te1(k), te2(k), s1, ds1
                   if (prnt2d) write (31, 2500) te1(k), te2(k), s1, ds1
                endif
             endif
             if (ipar < 7) then
                s5(k) = d2spesp(i4,i,k)*dtdo(k,i)
                ds5(k) = sqrt(abs(d2spesp(i4,i,k)))*dtdo(k,i)
                s6(k) = d2spepf(i4,i,k)*dtdo(k,i)
                ds6(k) = sqrt(abs(d2spepf(i4,i,k)))*dtdo(k,i)
                s7(k) = d2spef(i4,i,k)*dtdo(k,i)
                ds7(k) = sqrt(abs(d2spef(i4,i,k)))*dtdo(k,i)
                ls1 = s1 > zro .or. s2 > zro .or. s3 > zro .or. &
                     & s4 > zro
                if (ls1) then
                   if (prnt1d) write (31, 2400) te1(k), te2(k), s1, ds1, &
                        & s2, ds2, s3, ds3, s4, ds4
                   if (prnt2d) write (31, 2500) te1(k), te2(k), s1, ds1, &
                        & s2, ds2, s3, ds3, s4, ds4
                endif
             endif
          end do
!   End of k loop (energy for double differential cross sections)
          s1 = d2spec(i4,i,nt3)*doo(i)
          ds1 = sqrt(abs(d2spec(i4,i,nt3)))*doo(i)
          if (ipar < 7) then
             if (ipar < 3) then
                s2 = d2spe(i1,i,nt3)*doo(i)
                ds2 = sqrt(abs(d2spe(i1,i,nt3)))*doo(i)
                s3 = d2spe(i2,i,nt3)*doo(i)
                ds3 = sqrt(abs(d2spe(i2,i,nt3)))*doo(i)
                s4 = d2spe(i3,i,nt3)*doo(i)
                ds4 = sqrt(abs(d2spe(i3,i,nt3)))*doo(i)
             elseif (ipar >= 3) then
                s2 = d2speco(i5,i,nt3)*doo(i)
                ds2 = sqrt(abs(d2speco(i5,i,nt3)))*doo(i)
                s3 = d2spe(i1,i,nt3)*doo(i)
                ds3 = sqrt(abs(d2spe(i1,i,nt3)))*doo(i)
                s4 = d2spe(i2,i,nt3)*doo(i)
                ds4 = sqrt(abs(d2spe(i2,i,nt3)))*doo(i)
             endif
          elseif (ipar >= 7) then
             write (31, 3000) s1, ds1
          endif
          if (ipar < 7) then
             ls1 = s1 > zro .or. s2 > zro .or. s3 > zro .or. &
                  & s4 > zro
             if (ls1) write (31, 3000) s1, ds1, s2, ds2, s3, ds3, s4, ds4
             if (output%doubleDiffSpectra == 2) then
                write (31, 1100) par(ipar)
                do k = 1,ntt
                   prnt1d = (dble(int(dtt(k)*10 + 1.d-6))/10.d0)  ==  dtt(k)
                   prnt2d = .not.prnt1d
                   ls2 = s5(k) > zro .or. s6(k) > zro .or. s7(k) > zro
                   if (ls2) then
                      if (prnt1d) write (31, 2400) te1(k), te2(k), s5(k), &
                           & ds5(k), s6(k), ds6(k), s7(k), ds7(k)
                      if (prnt2d) write (31, 2500) te1(k), te2(k), s5(k), &
                           & ds5(k), s6(k), ds6(k), s7(k), ds7(k)
                   endif
                end do
                s5a = d2spesp(i4,i,nt3)*doo(i)
                ds5a = sqrt(abs(d2spesp(i4,i,nt3)))*doo(i)
                s6a = d2spepf(i4,i,nt3)*doo(i)
                ds6a = sqrt(abs(d2spepf(i4,i,nt3)))*doo(i)
                s7a = d2spef(i4,i,nt3)*doo(i)
                ds7a = sqrt(abs(d2spef(i4,i,nt3)))*doo(i)
                ls2 = s5a > zro .or. s6a > zro .or. s7a > zro
                if (ls2) write (31, 3000) s5a, ds5a, s6a, ds6a, s7a, ds7a
             endif
          endif
       end do
!   End of i loop (theta for double differential cross sections)
    endif
    return

! ======================================================================
1000 format (/4x,'T',a,'[MeV]',12x,'Total',16x,'Cascade',15x, &
          & 'Precompound',8x,'Total Evaporation')
1100 format (/8x,'Components of Evaporation:'/ &
          & /4x,'T',a,'[MeV]',4x,'From evap. residues',9x, &
          & 'Prefission',10x,'Fission Fragments')
1300 format (/2x,'Ang.',a,'[deg.]',5x,'Total = Cascade')
1400 format (/4x,'T',a,'[MeV]',6x,'Total = Cascade')
1500 format (/2x,'Ang.',a,9x,'Total',16x,'Cascade',15x,'Precompound', &
          &  9x,'Total Evaporation'/2x,'[deg.]')
1550 format (/2x,'Ang.',a,9x,'Total',14x,'Coalescence',13x, &
          &  'Precompound',9x,'Total Evaporation'/2x,'[deg.]')
1600 format (/8x,'Components of Evaporation:'/ &
          & /2x,'Ang.',a,3x,'From evap. residues',8x,'Prefission', &
          & 10x,'Fission Fragments'/2x,'[deg.]')
1700 format (/4x,'T',a,'[MeV]',11x,'Total',15x,'Coalescence',12x, &
          & 'Precompound',9x,'Total Evaporation')
1800 format (/1x,'Double differential cross sections [mb/MeV/sr];'/1x, &
          & 'Lab. angle = ',f5.1,' to ',f5.1,' degrees.')
1900 format (/1x,'--------------------------- Energy Spectrum', &
          & ' [mb/MeV] --------------------------')
2000 format (/1x,'----------------Normalized Energy Probabili', &
          & 'ty Spectrum [1/MeV] ----------------')
2100 format (5x,'Energy spectrum from ',f5.1,' to ',f6.1,' MeV (zero ', &
          & 'values suppressed).')
2200 format (5x,'Energy spectrum from ',f6.2,' to ',f7.2,' MeV (zero ', &
          & 'values suppressed).')
2300 format (/1x,'------------------------ Angular Distributions ', &
          & '[mb/sr] ------------------------')
2400 format (1x,f6.1,'-',f6.1,1x,4(1x,1pe9.3,' +/- ',e8.2))
2500 format (1x,f6.2,'-',f6.2,1x,4(1x,1pe9.3,' +/- ',e8.2))
2700 format (2x,f5.1,2x,4(1x,1pe9.3,' +/- ',e8.2))
3000 format (/2x,'Integrated:  ',4(1x,1pe9.3,' +/- ',e8.2))
3200 format (/1x,'Integ.  ',4(1x,1pe9.3,' +/- ',e8.2))
5000 format ("Divide by zero error prevented in 'prtdist.f90', line(s) ", A)
! ======================================================================
  end subroutine prtdist
