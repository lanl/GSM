
   function tcul (fbuObj, k, ntv)

! ======================================================================
!
!    Calculates Coulomb energy (see Eqn. 73 in the CEM Manual)
!
!    Called by: RAZVAL
!
!     Last change: 13-Aug-2003 BY NVMokhov
!     Edited by A. J. Sierk, LANL T-16, September, 2003.
!     Edited by A. J. Sierk, LANL T-2, February, 2009.
!     Edited by AJS, LANL T-2, December, 2011.
!     Edited by LMK, XCP-3, July 2013 (included error protection).
!     Edited by CMJ, XCP-3, July 2018 (creation of FermiBreakup class)
!
!    k = Number of products
!    tcul = Coulomb barrier
!    ntv = Particle code (100*A + Z)
!    radncl = 1.5
!
! ======================================================================

     use, intrinsic :: iso_fortran_env, only: int32, real64
     use fermiBreakupParams, only : zro, one, thrd, twthrd, thousandth, &
          & ato3rd

    implicit none
     class(FermiBreakup), intent(inout) :: fbuObj
     integer(int32),      intent(in   ) :: k    ! Number of products
     integer(int32),      intent(in   ) :: ntv(fbuObj%options%recNumNucleons)   ! Particle code (100*A + Z)
     real(real64)                       :: tcul

     integer(int32) :: i, iaa
     real(real64)   :: a, ec, z, zi
     integer(int32), dimension(fbuObj%options%recNumNucleons) :: ia, iz

! ======================================================================

     real(real64), parameter :: r0     = 1.3_real64    ! Radius parameter
     real(real64), parameter :: kappa  = 1.0_real64    ! V/V_0 term
     real(real64), parameter :: eCoul2 = 1.44_real64   ! "e" variable in coulomb potential equation
     real(real64), parameter :: coef   = &
          & (0.6d0*eCoul2/r0)*(one/((one + kappa)**thrd))

! ======================================================================

     tcul = zro
     a = zro
     z = zro
     ! Obtain A, Z of each fragment in the break-up channel
     do i = 1,min(k, fbuObj%options%recNumNucleons)
        ia(i) = ntv(i)/100
        iz(i) = ntv(i) - ia(i)*100
        a = a + dble(ia(i))
        z = z + dble(iz(i))
     end do
     iaa = nint(a)

     ! Determine coulomb potential for each fragment (summation term)
     ec = zro
     do i = 1,k
        zi = dble(iz(i))
        ec = ec + zi*zi/ato3rd(ia(i))
     end do

     ! Subtract from coloumb potential of original nucleus the fragment's combined couloumb potential
     ec = z*z/ato3rd(iaa) - ec

     ! Obtain coulomb potential, ensure valid
     tcul = thousandth*( coef * ec )
     tcul = max (zro, tcul)

     return
! ======================================================================
   end function tcul
