
  subroutine prtmult (gsmObj, output, outData, ipar, fn, sigin, nt2, nt3)

! ======================================================================
!
!   This subroutine prints out multiplicities and average energies of
!   emitted particles of type ipar, including a breakdown into cascade,
!   preequilibrium, those evaporated from non fissioning systems, those
!   evaporated prior to fission, from fission fragments, and those
!   from coalescence of cascade nucleons.
!
!   Extracted from old TYPEOUT by A. J. Sierk, LANL T-16, December,
!   2003.
!   Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use OutputDataMod, only: OutputData
    use gsm_params, only: zro

    implicit none
    class(GSM),        intent(inout) :: gsmObj
    class(GSMOutput),  intent(in   ) :: output
    class(OutputData), intent(in   ) :: outData
    integer(int32),    intent(in   ) :: ipar
    real(real64),      intent(in   ) :: fn
    real(real64),      intent(in   ) :: sigin
    integer(int32),    intent(in   ) :: nt2
    integer(int32),    intent(in   ) :: nt3

    integer(int32) :: i, i1, i2, i3, i4, i5, i6, ipart, j
    real(real64)   :: tempfn
    character(len=5) :: partx = ""

! ======================================================================

    integer(int32), parameter, dimension(9, 5) :: ipp = reshape( &
         & ([  1,  2,  3,  1,  1, &
         &     4,  5,  6,  2,  5, &
         &     7,  7,  8,  3,  9, &
         &     9,  9, 10,  4, 12, &
         &    11, 11, 12,  5, 15, &
         &    13, 13, 14,  6, 18, &
         &     7,  7,  7,  7, 21, &
         &     8,  8,  8,  8, 22, &
         &     9,  9,  9,  9, 23  ] ), shape(ipp), order=([2, 1])  )

    character(LEN=5), parameter, dimension(23) :: part = &
         [ 'T  n ', 'C  n ', 'P  n ', 'E  n ', 'T  p ', 'C  p ', &
         & 'P  p ', 'E  p ', 'T  d ', 'P  d ', 'E  d ', 'T  t ', &
         & 'P  t ', 'E  t ', 'T He3', 'P He3', 'E He3', 'T He4', &
         & 'P He4', 'E He4', 'pi-  ', 'pi0  ', 'pi+  '           ]
    character(LEN=5), parameter, dimension( 6) :: partf = &
         & ['F  n ', 'F  p ', 'F  d ', 'F  t ', 'F he3', 'F he4']
    character(LEN=5), parameter, dimension( 6) :: partpf = &
         & ['Pf n ', 'Pf p ', 'Pf d ', 'Pf t ', 'Pfhe3', 'Pfhe4']
    character(LEN=5), parameter, dimension( 6) :: partsp = &
         & ['Sp n ', 'Sp p ', 'Sp d ', 'Sp t ', 'Sphe3', 'Sphe4']
    character(LEN=5), parameter, dimension( 4) :: pncoa = &
         & ['Co d ', 'Co t ', 'CoHe3', 'CoHe4']

! ======================================================================

    real(real64)   :: pmul, dpmul, y, dy, emean, pmulc, dpmulc, yc, &
         & dyc, emeanc, pmulpf, dpmulpf, ypf, dypf, emeanpf, &
         & pmulsp, dpmulsp, ysp, dysp, emeansp, pmulf, dpmulf, yf, dyf, &
         & emeanf, pmcoa, dpmcoa, ycoal, dycoal, emeanco
    common /multip/ pmul(9),   dpmul(9),   y(9),   dy(9),   emean(9), &
         & pmulc(14), dpmulc(14), yc(14), dyc(14), emeanc(14), &
         & pmulpf(6), dpmulpf(6), ypf(6), dypf(6), emeanpf(6), &
         & pmulsp(6), dpmulsp(6), ysp(6), dysp(6), emeansp(6), &
         & pmulf(6),  dpmulf(6),   yf(6),  dyf(6),  emeanf(6), &
         & pmcoa(4), dpmcoa(4), ycoal(4), dycoal(4), emeanco(4)
    real(real64)   :: speco, angco, d2speco
    common /rescoa/ speco(4,200), angco(4,20), d2speco(4,10,200)
    real(real64)   :: spef, angf, d2spef
    common /resfis/ spef(6,200), angf(6,20), d2spef(6,10,200)
    real(real64)   :: spe, an
    common /resmac/ spe(14,200), an(14,20)
    real(real64)   :: specsp, angsp, d2spesp
    common /respal/ specsp(6,200), angsp(6,20), d2spesp(6,10,200)
    real(real64)   :: specpf, angpf, d2spepf
    common /resprf/ specpf(6,200), angpf(6,20), d2spepf(6,10,200)
    real(real64)   :: spec, ang, chan, dadz
    common /result/ spec(9,200), ang(9,20), chan(218), dadz(351,151)

! ======================================================================

    i1 = ipp(ipar,1)
    i2 = ipp(ipar,2)
    i3 = ipp(ipar,3)
    i4 = ipp(ipar,4)
    i5 = ipp(ipar,5)
    i6 = ipar - 2

    tempfn = fn
    if ( tempfn < div0Lim .and. tempfn > -div0Lim ) then
       tempfn = div0Lim
       write(gsmObj%io%message,5000) "82"
       call gsmObj%io%print(4, 3, gsmObj%io%message)
    end if

    if (ipar == 1) then
       do i = 1,14
          if (spe(i,nt3) > zro) then
             emeanc(i) = spe(i,nt2)/spe(i,nt3)
             pmulc(i) = spe(i,nt3)/tempfn
             dpmulc(i) = sqrt(abs(spe(i,nt3)))/tempfn
             yc(i) = pmulc(i)*sigin
             dyc(i) = dpmulc(i)*sigin
          else
             emeanc(i) = zro
             pmulc(i) = zro
             dpmulc(i) = zro
             yc(i) = zro
             dyc(i) = zro
          endif
       end do
       do j = 1,9
          if (spec(j,nt3) > zro) then
             emean(j) = spec(j,nt2)/spec(j,nt3)
             pmul(j) = spec(j,nt3)/tempfn
             dpmul(j) = sqrt(abs(spec(j,nt3)))/tempfn
             y(j) = pmul(j)*sigin
             dy(j) = dpmul(j)*sigin
          else
             emean(j) = zro
             pmul(j) = zro
             dpmul(j) = zro
             y(j) = zro
             dy(j) = zro
          endif
       end do
       do i = 1,4
          pmcoa(i)  = outData%pcoal(i)/tempfn
          dpmcoa(i) = sqrt(abs(outData%pcoal(i)))/tempfn
          ycoal(i)  = pmcoa(i)*sigin
          dycoal(i) = dpmcoa(i)*sigin
          if (speco(i,nt3) > zro) then
             emeanco(i) = speco(i,nt2)/speco(i,nt3)
          else
             emeanco(i) = zro
          endif
       end do
       write (31, 1000) outData%ncoal
       do  i = 1,6
          pmulf(i) = spef(i,nt3)/tempfn
          dpmulf(i) = sqrt(abs(spef(i,nt3)))/tempfn
          yf(i) = pmulf(i)*sigin
          dyf(i) = dpmulf(i)*sigin
          if (spef(i,nt3) > zro) then
             emeanf(i) = spef(i,nt2)/spef(i,nt3)
          else
             emeanf(i) = zro
          endif
          pmulsp(i) = specsp(i,nt3)/tempfn
          dpmulsp(i) = sqrt(abs(specsp(i,nt3)))/tempfn
          ysp(i) = pmulsp(i)*sigin
          dysp(i) = dpmulsp(i)*sigin
          if (specsp(i,nt3) > zro) then
             emeansp(i) = specsp(i,nt2)/specsp(i,nt3)
          else
             emeansp(i) = zro
          endif
          pmulpf(i) = specpf(i,nt3)/tempfn
          dpmulpf(i) = sqrt(abs(specpf(i,nt3)))/tempfn
          ypf(i) = pmulpf(i)*sigin
          dypf(i) = dpmulpf(i)*sigin
          if (specpf(i,nt3) > zro) then
             emeanpf(i) = specpf(i,nt2)/specpf(i,nt3)
          else
             emeanpf(i) = zro
          endif
       end do
       write (31, 1100)
       if (output%fisOnly) write (31, 1200)
       write (31, 1300)
       write (31, 1400)
    endif
    if (pmul(ipar).ne.zro .and. ipar < 7) then
       if (ipar > 1) write (31, 1400)
       ipart = i5
       partx = part(ipart)
!  Total:
       write (31, 1500) partx, pmul(i4), dpmul(i4), y(i4), dy(i4), &
            & emean(i4)
       ipart = ipart + 1
       partx = part(ipart)
!  Cascade:
       if (pmulc(i1).ne.zro .and. ipar < 3) &
            & write (31, 1500) partx, pmulc(i1), dpmulc(i1), yc(i1), &
            & dyc(i1), emeanc(i1)
       if (ipar < 3) ipart = ipart + 1
       partx = part(ipart)
!  Preequilibrium:
       if (pmulc(i2).ne.zro) &
            & write (31, 1500) partx, pmulc(i2), dpmulc(i2), yc(i2), &
            & dyc(i2), emeanc(i2)
       ipart = ipart + 1
       partx = part(ipart)
!  Evaporation from non-fissioning systems:
       if (pmulsp(i4).ne.zro) &
            & write (31, 1500) partsp(i4), pmulsp(i4), dpmulsp(i4), ysp(i4), &
            & dysp(i4), emeansp(i4)
!  Evaporation from fissioning systems:
       if (pmulpf(i4).ne.zro) &
            & write (31, 1500) partpf(i4), pmulpf(i4), dpmulpf(i4), ypf(i4), &
            & dypf(i4), emeanpf(i4)
!  Evaporation from fission fragments:
       if (pmulf(i4).ne.zro) &
            & write (31, 1500) partf(i4), pmulf(i4), dpmulf(i4), yf(i4), &
            & dyf(i4), emeanf(i4)
!  Total Evaporation:
       if (pmulc(i3).ne.zro) &
            & write (31, 1500) partx, pmulc(i3), dpmulc(i3), yc(i3), &
            & dyc(i3), emeanc(i3)
!  Coalescence:
       if (ipar > 2 .and.i6 >= 1)  then
          if(pmcoa(i6).ne.zro) &
               & write (31, 1500) pncoa(i6), pmcoa(i6), dpmcoa(i6), ycoal(i6), &
               & dycoal(i6), emeanco(i6)
       end if
    elseif (ipar >= 7) then
!  Pions:
       if (pmul(7).ne.zro .or. pmul(8).ne.zro .or. pmul(9).ne.zro) then
          if (ipar == 7) write (31, 1400)
          if (pmul(ipar).ne.zro) then
             if (pmul(i4).ne.zro) &
                  & write (31, 1500) part(i5), pmul(i4), dpmul(i4), y(i4), &
                  & dy(i4), emean(i4)
          endif
       endif
    endif
    if (ipar == 9) write (31, 1400)
    return

! ======================================================================
1000 format (/1x,'Number of coalesced d, t, He3, He4, He6, Li6, Li7, ', &
          & 'Be7 = ', 8i8)
1100 format (/1x,'Mean multiplicities, yields, and mean energies ', &
          & 'of ejected particles:'/1x,'(Notation: T - all production', &
          & ' mechanisms, C - cascade, P - pre-equilibrium,'/1x, &
          & 'Sp - from spallation residues, Pf - from nuclei before ', &
          & 'fission,'/1x,'F - from fission fragments, E - total ', &
          & 'evaporation = Sp + Pf + F,'/1x,'Co - Coalescence from ', &
          & 'cascade;'/1x,'Values which are identically zero are ', &
          & 'not printed.')
1200 format (/1x,'The multiplicities are printed for fission events ', &
          & 'only!')
1300 format (/1x,'Part.',5x,'Multiplicities',11x,'Yields [mb]', &
          & 5x,'<TKE> [MeV]')
1400 format (1x,'****************************************************', &
          & '*************')
1500 format (1x,a,f10.4,' +/- ',f6.4,2x,f10.3,' +/- ',f7.3,2x,f8.2)
5000 format ("Divide by zero error prevented in 'prtmult.f90', line(s) ", A)
! ======================================================================
  end subroutine prtmult
