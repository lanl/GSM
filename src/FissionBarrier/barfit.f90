
  subroutine barfit (fbObj, a, z, il, bfis, egs)

! ======================================================================
!
!***********************************************************************
!          S. G. Mashnik intoduced (25.10.1992) the following changes: *
!******************   z=zsierk, a=asierk, iz=z, ia=a   *****************
!                                                                      *
!   Removed unnecessary real*8 variables;                              *
!   Other changes made to fit into CEM95; A. J. Sierk, February, 1996. *
!   Made consistently real*8 for all f.p. variables;
!   Removed EQUIVALENCE statements;
!   Gives egs = egs(Lmax) for L >= Lmax.
!   Edited by A. J. Sierk, LANL T-16, October, 2003.
!   Edited by LMK, LANL XCP-3, July 2013 (included error protection).
!
!***********************************************************************
!
!    This subroutine returns the barrier height bfis, and the rotating
!    ground state energy, egs (in MeV) when called with real arguments
!    z, the atomic number, a, the atomic mass number, and (integer) il,
!    the angular momentum in units of h-bar, (Planck's constant divided
!    by 2*pi).
!
!         The fission barrier for il = 0 is calculated from a 7th order
!    fit in two variables to 638 calculated fission barriers for z 
!    values from 20 to 110.  These 638 barriers are fit with an rms 
!    deviation of 0.10 Mev by this 49-parameter function.
!    If BARFIT is called with (z,a) values outside the range of the
!    fit the barrier height is set to 0.0, and a message is
!    printed on unit 16.
!
!         For il values not equal to zero, the values of
!    l at which the barrier is  80%  and  20%  of the l=0 value are
!    respectively fit to 20-parameter functions of  z  and  a, over a
!    more restricted range of  a  values, than is the case for  l = 0.
!    The value of l where the barrier disappears, lmax, for 81 nuclei,
!    is fit to a 40-parameter function of z and a.
!         Once again, if a (z,a) pair is outside of the range of
!    validity of the fit, the barrier value is set to 0.0 and a message
!    is printed. These three values  (bfis(l=0), l-80, and l-20) and the
!    constraints of bfis = 0 and  d(bfis)/dl = 0 at l = lmax and 
!    d(bfis)/dl = 0 at l = 0 lead to a fifth-order fit to bfis(l) for
!    l>= l-20.  The last constraint, and the values l-80 and l-20,
!    and bfis(l=0) lead to a third-order fit for the region l <= l-20.
!
!         The ground-state energies are calculated from a 175-parameter
!    fit in Z, A, and L to 329 ground-state energies for 36 different
!    Z  and  A  values.
!    (The range of Z and A is the same as for L-80, L-20, and L-max)
!
!         The calculated barriers from which the fits were
!    made were calculated in 1983-1985 by A. J. Sierk of Los Alamos
!    National Laboratory, Group T-9, using  Yukawa-plus-exponential 
!    double folded nuclear energy, exact Couloub diffuseness corrections,
!    and diffuse-matter moments of inertia. The parameters of the model
!    are those derived by Moller and Nix in 1979:
!    r-0 = 1.16 fm, as = 21.13 Mev, kappa-s = 2.3  a = 0.68 fm.
!    The diffuseness of the matter and charge distributions used
!    corresponds to a surface diffuseness parameter (defined by Myers)
!    of 0.99 fm.  The fitted barriers for l = 0 are deviate from the 
!    calculated ones by less than 0.1 MeV;  The output from this
!    subroutine is a little less accurate.  Worst errors may be as large
!    as 0.5 MeV; characteristic uncertainty is in the range of 0.1-0.2
!    MeV.  The values of egs are generally approximated to within
!    about 0.1-0.2 MeV;  the largest deviation is about 0.5 Mev,
!    near L-I for light nuclei.

!         The rms deviation of lmax from the 81 input values is ~0.1
!    h-bar.  The approximate value found by this code is nearly always
!    within 0.5 h-bar of the calculated one.
!
!    Below is a table of test values to check implementation
!    of the program:
!               z, a, il  fiss bar lmax    egs
!              
!              28, 58, 0   33.14   45.8   0.00
!                    ,25   19.50         21.61
!                    ,40    2.88         50.12
!                    ,45.   0.07         58.18
!              65,153, 0   28.88   82.1   0.00
!                    ,50   16.16         19.07
!                    ,80    0.23         45.37
!                    ,82.   0.00         46.95
!              93,229, 0    3.76   68.1   0.00
!                    ,45    1.26          8.22
!                    ,68.   0.00         17.93
!
! ======================================================================
!
!    written by A. J. Sierk,  LANL  T-9
!    version 1.0   February, 1984
!    version 1.1   January, 1985:  improved coefficients in egs and lmax
!    version 1.2   September, 1985:  improved lmax, egs coefficients
!    version 1.21  June, 1986:   minor changes made
!    version 1.3   February, 1996: Moved elmax to subroutine, improved
!                  elmax coefficients.
!
!    Modified to calculate approximate egs for L > Lmax, A > 250, or
!       Z > 100;  AJS  02/14/05.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!
! ======================================================================
!
!      Copyright, 1996,  the Regents of the University of California.
!      This software was produced under a U. S. government contract
!      (w-7405-eng-36) by the Los Alamos National Laboratory, which is
!      operated by the University of California for the U. S. Department
!      of Energy.  The U. S. government is licensed to use, reproduce,
!      and distribute this software.  Permission is granted to the public
!      to copy and use this software without charge, provided that this
!      notice and any statement of authorship are reproduced on all
!      copies.  Neither the government nor the University makes any
!      warranty, expressed or implied, or assumes any liability
!      or responsibility for the use of this software.
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use fissionBarrierParams, only: zro, one, two, thr, four

    implicit none
    class(FissionBarrier), intent(inout) :: fbObj
    real(real64),          intent(in   ) :: a
    real(real64),          intent(in   ) :: z
    integer(int32),        intent(in   ) :: il
    real(real64),          intent(  out) :: egs
    real(real64),          intent(  out) :: bfis

    integer(int32) :: i, iz, j, k, l, m
    real(real64)   :: a1, a2, aa, aj, ak, amax, amin, bfdum, el, el20, &
         & el80, ell, elmax, q, qa, qb, temp, x, y, z1, zz

    real(real64), dimension( 7) :: pa = zro, pz= zro
    real(real64), dimension(10) :: pl = zro

! ======================================================================

    real(real64), dimension(5,4) :: emncof
    data emncof  & 
         /-9.01100d+2,-1.40818d+3, 2.77000d+3,-7.06695d+2, 8.89867d+2, &  
         1.35355d+4,-2.03847d+4, 1.09384d+4,-4.86297d+3,-6.18603d+2, & 
         -3.26367d+3, 1.62447d+3, 1.36856d+3, 1.31731d+3, 1.53372d+2, &  
         7.48863d+3,-1.21581d+4, 5.50281d+3,-1.33630d+3, 5.05367d-2/

    real(real64), dimension(5,4) :: elmcof
    data elmcof  & 
         /1.84542d+3,-5.64002d+3, 5.66730d+3,-3.15150d+3, 9.54160d+2, & 
         -2.24577d+3, 8.56133d+3,-9.67348d+3, 5.81744d+3,-1.86997d+3, &  
         2.79772d+3,-8.73073d+3, 9.19706d+3,-4.91900d+3, 1.37283d+3, & 
         -3.01866d+1, 1.41161d+3,-2.85919d+3, 2.13016d+3,-6.49072d+2/

    real(real64), dimension(7,7) :: elzcof
    data elzcof  &
         /5.11819909d+5,-1.30303186d+6, 1.90119870d+6,-1.20628242d+6, &  
         5.68208488d+5, 5.48346483d+4,-2.45883052d+4, & 
         -1.13269453d+6, 2.97764590d+6,-4.54326326d+6, 3.00464870d+6, & 
         -1.44989274d+6,-1.02026610d+5, 6.27959815d+4, &  
         1.37543304d+6,-3.65808988d+6, 5.47798999d+6,-3.78109283d+6, &  
         1.84131765d+6, 1.53669695d+4,-6.96817834d+4, & 
         -8.56559835d+5, 2.48872266d+6,-4.07349128d+6, 3.12835899d+6, & 
         -1.62394090d+6, 1.19797378d+5, 4.25737058d+4, &  
         3.28723311d+5,-1.09892175d+6, 2.03997269d+6,-1.77185718d+6, &  
         9.96051545d+5,-1.53305699d+5,-1.12982954d+4, &  
         4.15850238d+4, 7.29653408d+4,-4.93776346d+5, 6.01254680d+5, & 
         -4.01308292d+5, 9.65968391d+4,-3.49596027d+3, & 
         -1.82751044d+5, 3.91386300d+5,-3.03639248d+5, 1.15782417d+5, & 
         -4.24399280d+3,-6.11477247d+3, 3.66982647d+2/

    real(real64), dimension(5,7,5) :: egscof
    data ((egscof(i,j,1),i=1,5),j=1,7) &  
         /-1.781665232d6,-2.849020290d6, 9.546305856d5, 2.453904278d5, 3.656148926d5, &    
         4.358113622d6, 6.960182192d6,-2.381941132d6,-6.262569370d5, -9.026606463d5, &   
         -4.804291019d6,-7.666333374d6, 2.699742775d6, 7.415602390d5, 1.006008724d6, &    
         3.505397297d6, 5.586825123d6,-2.024820713d6,-5.818008462d5, -7.353683218d5, &   
         -1.740990985d6,-2.759325148d6, 1.036253535d6, 3.035749715d5, 3.606919356d5, &    
         5.492532874d5, 8.598827288d5,-3.399809581d5,-9.852362945d4, -1.108872347d5, &   
         -9.229576432d4,-1.431344258d5, 5.896521547d4, 1.772385043d4, 1.845424227d4/

    data ((egscof(i,j,2),i=1,5),j=1,7) & 
         /4.679351387d6, 7.707630513d6,-2.718115276d6,-9.845252314d5, -1.107173456d6, & 
         -1.137635233d7,-1.870617878d7, 6.669154225d6, 2.413451470d6, 2.691480439d6, &  
         1.237627138d7, 2.030222826d7,-7.334289876d6,-2.656357635d6,-2.912593917d6, & 
         -8.854155353d6,-1.446966194d7, 5.295832834d6, 1.909275233d6, 2.048899787d6, &  
         4.290642787d6, 6.951223648d6,-2.601557110d6,-9.129731614d5,-9.627344865d5, & 
         -1.314924218d6,-2.095971932d6, 8.193066795d5, 2.716279969d5, 2.823297853d5, &  
         2.131536582d5, 3.342907992d5,-1.365390745d5,-4.417841315d4, -4.427025540d4/

    data ((egscof(i,j,3),i=1,5),j=1,7) &
         /-3.600471364d6,-5.805932202d6, 1.773029253d6, 4.064280430d5, 7.419581557d5, &  
         8.829126250d6, 1.422377198d7,-4.473342834d6,-1.073350611d6, -1.845960521d6, & 
         -9.781712604d6,-1.575666314d7, 5.161226883d6, 1.341287330d6, 2.083994843d6, &  
         7.182555931d6, 1.156915972d7,-3.941330542d6,-1.108259560d6, -1.543982755d6, & 
         -3.579820035d6,-5.740079339d6, 2.041827680d6, 5.981648181d5, 7.629263278d5, &  
         1.122573403d6, 1.777161418d6,-6.714631146d5,-1.952833263d5, -2.328129775d5, & 
         -1.839672155d5,-2.871137706d5, 1.153532734d5, 3.423868607d4, 3.738902942d4/

    data ((egscof(i,j,4),i=1,5),j=1,7) &
         / 2.421750735d6, 4.107929841d6,-1.302310290d6,-5.267906237d5, -6.197966854d5, & 
         -5.883394376d6,-9.964568970d6, 3.198405768d6, 1.293156541d6, 1.506909314d6, &  
         6.387411818d6, 1.079547152d7,-3.517981421d6,-1.424705631d6,-1.629099740d6, & 
         -4.550695232d6,-7.665548805d6, 2.530844204d6, 1.021187317d6, 1.141553709d6, &  
         2.182540324d6, 3.646532772d6,-1.228378318d6,-4.813626449d5, -5.299974544d5, & 
         -6.518758807d5,-1.070414288d6, 3.772592079d5, 1.372024952d5, 1.505359294d5, &  
         9.952777968d4, 1.594230613d5,-6.029082719d4,-2.023689807d4,-2.176008230d4/

    data ((egscof(i,j,5),i=1,5),j=1,7) &
         /-4.902668827d5,-8.089034293d5, 1.282510910d5,-1.704435174d4, 8.876109934d4, &  
         1.231673941d6, 2.035989814d6,-3.727491110d5, 4.071377327d3,-2.375344759d5, & 
         -1.429330809d6,-2.376692769d6, 5.216954243d5, 7.268703575d4, 3.008350125d5, &  
         1.114306796d6, 1.868800148d6,-4.718718351d5,-1.215904582d5, -2.510379590d5, & 
         -5.873353309d5,-9.903614817d5, 2.742543392d5, 9.055579135d4, 1.364869036d5, &  
         1.895325584d5, 3.184776808d5,-9.500485442d4,-3.406036086d4, -4.380685984d4, & 
         -2.969272274d4,-4.916872669d4, 1.596305804d4, 5.741228836d3, 6.669912421d3/

! ======================================================================
!
!    The program starts here
! 
! ======================================================================

    iz = nint(z)
    el = dble(il)
    if (iz < 19 .or. iz > 111) then
       bfis = fbObj%bf2 (a, z, il, egs)
       return
    endif
    if (iz > 102 .and. il > 0) then
       bfis = fbObj%bf2 (a, z, il, egs)
       return
    endif
    amin = 1.2d0*z + 0.01d0*z*z
    amax = 5.8d0*z - 0.024d0*z*z
    if (a < amin .or. a > amax) then
       bfis = fbObj%bf2 (a, z, il, egs)
       return
    endif
    aa = 2.5d-3*a
    zz = 1.d-2*z
    bfis = zro
    call lpoly2 (zz, pz, 7)
    call lpoly2 (aa, pa, 7)
    do i = 1,7
       do j = 1,7
          bfis = bfis + elzcof(j,i)*pz(j)*pa(i)
       end do
    end do
    egs = zro
    if (il < 1) return
    amin = 1.4d0*z + 0.009d0*z*z
    amax = 20.d0 + thr*z
    if (a < amin-5.d0 .or. a > amax+10.d0) then
       bfis = fbObj%bf2 (a, z, il, egs)
       return
    endif
    el80 = zro
    el20 = zro
    do i = 1,4
       do j = 1,5
          el80 = el80 + elmcof(j,i)*pz(j)*pa(i)
          el20 = el20 + emncof(j,i)*pz(j)*pa(i)
       end do
    end do
    call elmaxc (z, a, elmax)
    if (el > elmax) then
       bfis = zro
       z1 = one
    else 
       temp = elmax
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(fbObj%io%message,1000) "272-274"
          call fbObj%io%print(4, 3, fbObj%io%message)
       end if
       x = el20/temp
       y = el80/temp
       z1 = el/temp
       if (el <= el20) then
          temp = el20**2*el80**2*(el20 - el80)
          if (temp < div0Lim .and. temp > -div0Lim) then
             temp = div0Lim
             write(fbObj%io%message,1000) "281"
             call fbObj%io%print(4, 3, fbObj%io%message)
          end if
          q = 0.20d0/(temp)
          qa =  q*(four*el80**3 - el20**3)
          qb = -q*(four*el80**2 - el20**2)
          bfis = bfis*(one + qa*el**2 + qb*el**3)
       else
          aj = (-20.d0*x**5 + 25.d0*x**4 - four)*(y - one)**2*y*y
          ak = (-20.d0*y**5 + 25.d0*y**4 - one)*(x - one)**2*x*x
          temp = (y - x)*((one - x)*(one - y)*x*y)**2
          if (temp < div0Lim .and. temp > -div0Lim) then
             temp = div0Lim
             write(fbObj%io%message,1000) "293"
             call fbObj%io%print(4, 3, fbObj%io%message)
          end if
          q = 0.2d0/(temp)
          qa =  q*(aj*y - ak*x)
          qb = -q*(aj*(two*y + one) - ak*(two*x + one))
          a1 = four*z1**5 - 5.d0*z1**4 + one
          a2 = qa*(two*z1 + one)
          bfis = bfis*(a1 + (z1 - one)*(a2 + qb*z1)*z1*z1*(z1 - one))
       endif
       bfis = max(bfis,zro)
    endif

! Now calculate rotating ground-state energy
    ell = z1
    call lpoly2 (ell, pl, 9)
    egs = zro
    if (el <= elmax .and. a <= 250.d0 .and. z <= 100.d0) then
       do k = 1,5
!   Use powers of (A/250) from 0-4.
          do l = 1,7
!   Use powers of (Z/100) from 0-6.
             do m = 1,5
!   Use even powers of (l/lmax) from 0-8.
                egs = egs + egscof(m,l,k)*pz(l)*pa(k)*pl(2*m-1)
             end do
          end do
       end do
    else
       bfdum = fbObj%bf2 (a, z, il, egs)
    endif
    egs = max(egs, zro)

    return
! ======================================================================
1000 format("Divide by zero error prevented in 'barfit.f90', line(s) ", A)
! ======================================================================
  end subroutine barfit


  subroutine elmaxc (zn, an, elmax)
!  subroutine elmaxc (fbObj, zn, an, elmax)

! ======================================================================
!   This subroutine evaluates a 2-dimensional fitting function
!   to arrive at an approximation to Lmax(Z,A), the value of
!   angular momentum in units of h-bar at which the fission
!   saddle point and the fission barrier vanishes.  Lmax is the
!   maximum amount of angular momentum which the nucleus can
!   sustain without fissioning.
!
!   The 2-dimensional fit was made to 192 calculated values of
!   Lmax for values of Z from 10 to 110, with A values spanning
!   a range larger than any nuclei which might be formed in
!   a reaction.
!
!   The 192 calculated values of Lmax were calculated in the
!   three-quadratic-surface axially symmetric parametrization,
!   or in the triaxial Legendre Polynomial parametrization
!   using Yukawa-plus-exponential nuclear energy, diffuse
!   Coulomb energy, and diffuse-matter moments of inertia by
!   A. J. Sierk of Los Alamos National Laboratory, Group T-2.
!
!   This subroutine written February, 1994.
!   Updated constants and extendd Z and A range of fit, May, 1996.
!   Edited by A. J. Sierk, LANL T-16, October, 2003.
!   Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64

    implicit none
!    class(FissionBarrier), intent(inout) :: fbObj   ! UNUSED!
    real(real64),          intent(in   ) :: zn
    real(real64),          intent(in   ) :: an
    real(real64),          intent(  out) :: elmax

    integer(int32) :: ia, ijb, iz
    real(real64)   :: xa, xz

    real(real64), dimension(npa) :: pla
    real(real64), dimension(npz) :: plz

! ======================================================================

    real(real64), parameter, dimension(40) :: b = &
         & [   -3.30787006d5,  8.34310803d5, -4.11955431d5,  9.82177374d5, &     
         &      3.40653902d5,  8.39479538d5, -2.13869250d6,  1.11540268d6, &    
         &     -2.52645986d6, -7.72016242d5, -1.00538896d6,  2.59234090d6, &    
         &     -1.49269622d6,  3.06066515d6,  7.02514732d5,  7.97067135d5, &    
         &     -2.08607738d6,  1.33444088d6, -2.52998737d6, -3.45973915d5, &    
         &     -4.15169782d5,  1.10638461d6, -7.71715934d5,  1.45998634d6, &     
         &      6.29116062d4,  1.07407080d5, -3.04281993d5,  2.18876998d5, &    
         &     -5.37873263d5,  2.14561453d4,  1.64345870d3,  6.52424621d3, &     
         &      3.41057615d3,  1.03029197d5, -1.17188782d4, -5.44933809d3, &     
         &      1.23028722d4, -1.25511319d4, -5.30918539d3,  1.17305352d3   ]

! ======================================================================

    xz = zn/100.d0
    xa = an/320.d0
    elmax = 0.d0

    ! Obtain legendre polynomials
    call lpoly2 (xz, plz, npz)
    call lpoly2 (xa, pla, npa)

    ! Determine fission barrier saddle point
    do iz = 1,npz
       do ia = 1,npa
          ijb = ia + npa*(iz-1)
          elmax = elmax + b(ijb)*plz(iz)*pla(ia)
       end do
    end do

    return
! ======================================================================
  end subroutine elmaxc


  subroutine lpoly2 (x, ple, n)
!  subroutine lpoly2 (fbObj, x, ple, n)

! ======================================================================
!
!  Calculates the first n Legendre polynomials, those of order 0 to n-1.
!
!   Written by A. J. Sierk,  LANL T-2, February, 1994.
!   Modified by A. J. Sierk, LANL T-16, October, 2003.
!   Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use fissionBarrierParams, only: hlf, one, two, thr, fiv, seven

    implicit none
!    class(FissionBarrier), intent(inout)               :: fbObj   ! UNUSED!
    integer(int32),        intent(in   )               :: n
    real(real64),          intent(in   )               :: x
    real(real64),          intent(out  ), dimension(n) :: ple

    integer(int32) :: i
    real(real64)   :: dim1, dim2, d2im3

! ======================================================================

    ! The first polynomials (1 through 4) are:
    if ( n >= 1 ) ple(1) = one
    if ( n >= 2 ) ple(2) = x
    if ( n >= 3 ) ple(3) = hlf*(thr*x*x - one)
    if ( n >= 4 ) ple(4) = hlf*x*(fiv*x*x - thr)

    dim2 = thr
    dim1 = two*two
    d2im3 = seven

    ! Get higher order polynomials (>= 5)
    if (n >= 5) then
       do i = 5,n
          ple(i) = (d2im3*x*ple(i-1) - dim2*ple(i-2))/dim1
          dim2 = dim1
          dim1 = dble(i)
          d2im3 = two + d2im3
       end do
    endif

    return
! ======================================================================
  end subroutine lpoly2


  subroutine bsfit (fbObj, zn, an, il, bs)

! ======================================================================
!
!   This subroutine evaluates a 3-dimensional fitting function
!   to arrive at an approximation to Bs(Z,A,L), the value of
!   the surface area of the macroscopic saddle-point shape in units
!   of the area of the spherical nucleus.
!
!   The 200-parameter 3-dimensional fit was made to 1053 calculated
!   values of Bs for values of Z from 15 to 100, with A values spanning
!   a range larger than any nuclei which might be formed in
!   a reaction, and for L from 0 to Lmax.
!
!   The 1053 calculated values of Bs were calculated in the
!   Legendre Polynomial axially symmetric parametrization,
!   or in the triaxial Legendre Polynomial parametrization
!   using Yukawa-plus-exponential nuclear energy, diffuse
!   Coulomb energy, and diffuse-matter moments of inertia by
!   A. J. Sierk of Los Alamos National Laboratory, Group T-2,
!   in May, 1996.
!
!   ELMAXC subroutine written February, 1994.
!   BSFIT modified from ELMAXC, May, 1996.
!   Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use fissionBarrierParams, only: zro

    implicit none
    class(FissionBarrier), intent(inout) :: fbObj
    real(real64),          intent(in   ) :: zn
    real(real64),          intent(in   ) :: an
    integer(int32),        intent(in   ) :: il
    real(real64),          intent(  out) :: bs

    integer(int32) :: ia, ijb, ill, iz, npl2
    real(real64)   :: el, elmax, temp, xa, xl, xz

    real(real64), dimension( npa)    :: pla = zro
    real(real64), dimension( npz)    :: plz = zro
    real(real64), dimension(9_int32) :: pll = zro

! ======================================================================

    real(real64), parameter, dimension(270) :: b = &
         & [  -17605.486d0,  7149.518d0, 23883.053d0, 21651.272d0,  -348.402d0, &
         &     67617.841d0,-43739.724d0,-72255.975d0,-79130.794d0,-44677.571d0, & 
         &    -38284.438d0,  1629.405d0, 49928.158d0, 41804.076d0, -4744.661d0, & 
         &      2124.997d0,-30258.577d0,-45416.684d0,-40347.682d0,-20718.220d0, &  
         &     -7009.141d0,-12274.997d0,  5703.797d0,  4375.711d0, -2519.378d0, &  
         &     44633.464d0,-19664.988d0,-59871.761d0,-55426.914d0,  -834.542d0, &
         &   -168222.448d0,111291.480d0,179939.092d0,198092.569d0,111347.893d0, &  
         &     97361.913d0,-9404.674d0,-15757.922d0,-108003.071d0,  7218.634d0, & 
         &    -79483.915d0, 77325.943d0,113433.935d0,100793.473d0, 51134.715d0, &  
         &     17974.287d0, 28084.202d0,-15187.024d0,-12244.758d0,  4662.638d0, & 
         &    -52775.944d0, 26372.246d0, 69061.418d0, 66412.802d0,  4456.150d0, & 
         &    191743.968d0,-131672.398d0,-204896.966d0,-227867.388d0, &
         &   -126609.907d0,-115704.145d0,  22053.324d0, 146315.778d0, & 
         &    131364.763d0,   1176.787d0,  89552.374d0, -92216.043d0, &
         &   -129901.381d0,-115369.145d0, -56853.433d0, -21565.916d0, & 
         &    -26595.626d0, 19490.177d0, 16758.541d0, -2102.079d0, 43089.492d0, & 
         &    -24205.277d0, -54145.651d0, -54916.885d0,  -6982.880d0, &
         &   -148893.417d0,105958.013d0,157835.761d0,178691.634d0, 97246.541d0, &  
         &     94988.919d0, -28100.912d0,-115999.016d0,-110379.031d0, & ! Items 1-89
         &    -10165.655d0, -68277.650d0, & ! Items 90-91
         &     75226.302d0,100957.904d0, 89635.423d0, 41887.718d0, 17720.511d0, &  
         &     14837.634d0,-17476.245d0,-15762.519d0, -1297.194d0,-25673.002d0, &  
         &     15584.030d0, 30291.575d0, 33008.500d0,  6161.347d0, 83671.721d0, & 
         &    -60540.161d0,-86493.251d0,-101211.404d0,-53525.729d0, & 
         &    -56872.290d0, 21837.845d0, 65632.326d0, 67209.922d0, 11367.245d0, &
         &     37417.011d0,-43832.496d0,-56023.110d0,-49995.410d0,-21521.767d0, & 
         &    -10406.752d0, -4483.654d0, 11314.537d0, 10435.311d0,  2206.541d0, &  
         &     11115.697d0, -6672.659d0,-11696.751d0,-14183.580d0, -3586.280d0, & 
         &    -33914.223d0, 23453.203d0, 32630.486d0, 40847.174d0, 21223.389d0, &  
         &     24704.794d0,-10226.568d0,-25497.313d0,-29082.833d0, -7190.865d0, & 
         &    -14734.982d0, 17457.735d0, 21457.946d0, 19706.266d0,  7745.219d0, &   
         &      4283.242d0,   214.103d0, -5058.390d0, -4715.987d0, -1303.557d0, &  
         &     -3196.271d0,  1579.967d0,  2718.777d0,  3997.275d0,  1361.541d0, &   
         &      9152.527d0, -5182.345d0, -7412.465d0,-10822.583d0, -5804.143d0, &  
         &     -7177.709d0,  2317.015d0,  5874.736d0,  8208.479d0,  2858.079d0, &   
         &      3919.329d0, -4068.562d0, -4982.296d0, -5065.564d0, -1948.636d0, &  
         &     -1147.488d0,   265.063d0,  1372.992d0,  1333.320d0,   416.369d0, &    
         &       442.331d0,  -219.051d0,  -337.663d0,  -575.255d0,  -232.634d0, & ! Items 92-180
         &     -1217.351d0,   590.169d0,   854.438d0,  1468.791d0,   811.481d0, &   
         &      1031.157d0,  -227.306d0,  -651.183d0, -1174.558d0,  -528.960d0, &   
         &      -536.534d0,   413.429d0,   533.977d0,   663.064d0,   278.203d0, &    
         &       154.548d0,   -54.623d0,  -169.817d0,  -183.006d0,   -63.999d0, &
         &         0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,   &
         &         0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,   &
         &         0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,   &
         &         0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,   &
         &         0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,   &
         &         0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,   &
         &         0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,   &
         &         0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,   &
         &         0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,   &
         &         0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,   &
         &         0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,   &
         &         0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,   &
         &         0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0,   &
         &         0.0d0,       0.0d0,       0.0d0,       0.0d0,       0.0d0   ] ! Items 181-270

! ======================================================================

    npl2 = 2*npl - 1
    xz = zn/100.d0
    xa = an/320.d0
    el = dble(il)
    call elmaxc (zn, an, elmax)
    temp = elmax
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(fbObj%io%message,1000) "569"
       call fbObj%io%print(4, 3, fbObj%io%message)
    end if
    xl = el/temp
    call lpoly2 (xz, plz, npz)
    call lpoly2 (xa, pla, npa)
    call lpoly2 (xl, pll, npl2)
    bs = zro
    do iz = 1,npz
       do ia = 1,npa
          do ill = 1,npl
             ijb = ill + npl*npa*(iz-1) + npl*(ia-1)
             bs = bs + b(ijb)*plz(iz)*pla(ia)*pll(2*ill-1)
          end do
       end do
    end do


    return
! ======================================================================
1000 format("Divide by zero error prevented in 'barfit.f90', line(s) ", A)
! ======================================================================
  end subroutine bsfit


  function bf2 (fbObj, a, z, ln, egs0)

! ======================================================================
!
!   Function to give fission barrier in Mev, given A, Z, L, and 
!   excitation energy at the saddle point, e. bf2 is for use in
!   BARFIT when the a Z, A, L values lie outside the range of the fit.
!   It gives an approximate value using the Krappe-Nix-Sierk model.
!   Results are fairly crude, but should not be needed too often,
!   except perhaps for Z < 19 and Z > 12 (just above Fermi breakup).
!
!   CEM95 written by S. G. Mashnik
!   Edited by A. J. Sierk,  LANL  T-2  February, 1996.
!   bf2 created from bf by AJS, December, 1997
!   Modified by AJS, March, 1999.
!   "Last" change: 12-AUG-2003 by NVMokhov
!   Modified by A. J. Sierk, LANL T-16, October, 2003.
!   Modified by AJS, February, 2005.
!   Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use fissionBarrierParams, only: zro, hlf, one, thr, emnucb, ato3rd

    implicit none
    class(FissionBarrier), intent(inout) :: fbObj
    real(real64),          intent(in   ) :: a
    real(real64),          intent(in   ) :: z
    integer(int32),        intent(in   ) :: ln
    real(real64),          intent(  out) :: egs0
    real(real64)                         :: bf2

    integer(int32) :: in, iz
    real(real64)   :: a2, a3, a3rt, a3rt2, acor, athrd, bf, bf0, bx, &
         & coulmb, daz, dbfi, emcorr, emscorr, emx, ergs0, ersad, fjrb, &
         & fjsad, gamma, sufnuc, sym, temp, um2, un, w, x, zsq

! ======================================================================

    un = a - z
    zsq = z*z
    temp = a
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(fbObj%io%message,1000) "642"
       call fbObj%io%print(4, 3, fbObj%io%message)
    end if
    sym = ((un - z)/temp)**2
    egs0 = 1.d-6
    iz = nint(z)
    in = nint(un)

    gamma = thr
    a2 = 21.7d0
    a3 = 0.7322d0
    athrd = ato3rd(nint(a))
    a3rt = athrd
    a3rt2 = a3rt**2
    acor = one - gamma*sym
    sufnuc = a2*acor*a3rt2
    coulmb = a3*zsq/a3rt
    temp = sufnuc
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(fbObj%io%message,1000) "666"
       call fbObj%io%print(4, 3, fbObj%io%message)
    end if
    x = hlf*coulmb/temp
    emx = fbObj%fbMolnix%defineEnergy (iz, in, 2)
    temp = emnucb*a
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(fbObj%io%message,1000) "673"
       call fbObj%io%print(4, 3, fbObj%io%message)
    end if
    emcorr = one + 0.001d0*emx/(temp)
    um2 = dble(ln)**2
!  Finite-range approximation to rigid sphere moment of inertia:
    fjrb = (fbObj%options%r0m*fbObj%options%r0m*athrd**2 + 10.d0*amm**2)*a
    temp = fjrb
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(fbObj%io%message,1000) "682"
       call fbObj%io%print(4, 3, fbObj%io%message)
    end if
    ergs0 = 52.252429d0*um2/temp
    temp = emcorr
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(fbObj%io%message,1000) "688"
       call fbObj%io%print(4, 3, fbObj%io%message)
    end if
    ergs0 = ergs0/temp
    if (x >= one) then
       bf2 = zro
    else
       bx = fbObj%subev (x, xx2, yy2, 51)
       bf0 = sufnuc*bx
       bf0 = max(bf0, zro)
       daz = fbObj%fbMolnix%shellEnergy(a, z)
       bf = bf0 - daz
!  Approximate rotational lowering of barrier:
!  w = Strutinsky approximation to LDM moment of inertia of saddle:
       w = fbObj%subev (x, xx1, ws, 51)
!  Finite-range approximation to rigid saddle moment of inertia:
!  (Diffuseness correction is independent of shape!!)
       fjsad = (w*fbObj%options%r0m*fbObj%options%r0m*athrd**2 + 10.d0*amm**2)*a
       temp = emnucb*a
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(fbObj%io%message,1000) "708"
          call fbObj%io%print(4, 3, fbObj%io%message)
       end if
       emscorr = one + 0.001d0*(emx + bf)/(temp)
       temp = fjsad
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(fbObj%io%message,1000) "714"
          call fbObj%io%print(4, 3, fbObj%io%message)
       end if
       ersad = 52.252429d0*um2/temp
       temp = emscorr
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(fbObj%io%message,1000) "720"
          call fbObj%io%print(4, 3, fbObj%io%message)
       end if
       ersad = ersad/temp
       dbfi = ersad - ergs0
       bf = bf + dbfi
       bf2 = max(bf, zro)
    endif
    egs0 = ergs0

    return
! ======================================================================
1000 format("Divide by zero error prevented in 'barfit.f90', line(s) ", A)
! ======================================================================
  end function bf2
