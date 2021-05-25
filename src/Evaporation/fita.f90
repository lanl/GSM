
  subroutine fitafac (atar0, ztar0, ener0, fitaf, fitaf1)

! ======================================================================
!
!      Quadratic interpolation of fitaf1 = f(z0, t0) for actinides.
!
!    Written by S. G. Mashnik
!   "Last" change: 12-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
!    Changed by K.K. Gudima, using new fit by M.I. Baznat, Nov., 2004
!    Edited by A. J. Sierk, LANL T-16, January, 2005.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64

    implicit none
    real(real64), intent(in   )  :: atar0
    real(real64), intent(in   )  :: ztar0
    real(real64), intent(in   )  :: ener0
    real(real64), intent(  out) :: fitaf
    real(real64), intent(  out) :: fitaf1

    integer(int32) :: ia, i, iz, k
    real(real64)   :: a, a1, b, b1, c, c1, d, da, da1, db, db1, dc, &
         & dc1, t0, x, x1, x12, x2, x23, x3, x31, y, y1, y2, &
         & y3, z, z0, z1, z2, z3

! ======================================================================

    integer(int32), parameter :: mq = 17

    ! For calculation
    real(real64),   dimension(mq) :: xq = 0.d0
    real(real64),   dimension(mq) :: yq = 0.d0
    real(real64),   dimension(mq) :: zq = 0.d0

    ! For fits
    real(real64),   dimension(5), parameter :: zfit = &
         & [ 89.d0, 90.d0, 92.d0, 93.d0, 94.d0 ]
    real(real64),   dimension(mq), parameter :: tfit = &
         & [   0.02d0, 0.03d0, 0.04d0, 0.05d0, 0.07d0, 0.1d0, 0.15d0, &
         &     0.2d0,  0.3d0,  0.4d0,  0.5d0,  0.8d0,  1.0d0, 1.5d0, &
         &     2.0d0,  2.5d0,  3.d0  ]

    ! Ac-227
    real(real64),   dimension(mq), parameter :: ac227 = &
         & [  1.20d0,  1.20d0,  1.20d0,  1.20d0,  1.20d0,  1.15d0,  1.022d0, &
         &    0.993d0, 0.965d0, 0.955d0, 0.944d0, 0.936d0, 0.946d0, 0.979d0, &
         &    0.995d0, 1.013d0, 1.18d0  ]
    real(real64),   dimension(mq), parameter :: ac227f1 = &
         & [   1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, &
         &     1.0d0, 1.0d0, 1.0d0, 1.0d0, 5.0d0, 1.0d0, 1.0d0, &
         &     1.0d0, 1.0d0, 1.0d0   ]

    ! Th-232
    real(real64),   dimension(mq), parameter :: th232 = &
         & [  1.0d0,   1.0d0,   1.0d0,   0.998d0, 0.997d0, 0.995d0, 0.994d0, &
         &    0.990d0, 0.978d0, 0.967d0, 0.950d0, 0.924d0, 0.933d0, 0.955d0, &
         &    0.971d0, 0.973d0, 0.974d0   ]
    real(real64),   dimension(mq), parameter :: th232f1 = &
         & [   1.35d0, 1.44d0, 1.44d0, 1.30d0, 1.25d0, 1.7d0, 1.9d0, &
         &     2.30d0, 3.70d0, 5.10d0, 6.90d0, 8.00d0, 9.7d0, 6.0d0, &
         &     4.00d0, 2.10d0, 1.00d0   ]

    ! U-233
    real(real64),   dimension(mq), parameter :: u233 = &
         & [   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0, &
         &     0.990d0, 0.981d0, 0.979d0, 0.973d0, 0.925d0, 0.931d0, 0.950d0, &
         &     0.975d0, 0.989d0, 1.007d0   ]
    real(real64),   dimension(mq), parameter :: u233f1 = &
         & [   4.30d0, 1.00d0, 1.00d0, 1.10d0, 1.40d0, 2.0d0, 2.1d0, &
         &     2.10d0, 3.50d0, 3.90d0, 6.50d0, 9.60d0,10.1d0, 4.2d0, &
         &     2.40d0, 1.50d0, 1.10d0   ]

    ! U-235
    real(real64),   dimension(mq), parameter :: u235 = &
         & [   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   0.998d0, 0.997d0, &
         &     0.997d0, 0.996d0, 0.993d0, 0.986d0, 0.950d0, 0.950d0, 0.966d0, &
         &     0.973d0, 0.991d0, 0.998d0   ]
    real(real64),   dimension(mq), parameter :: u235f1 = &
         & [   2.8d0, 1.0d0, 1.0d0, 1.0d0, 1.1d0, 1.6d0, 1.6d0, &
         &     1.7d0, 1.9d0, 3.0d0, 5.0d0, 9.6d0,10.2d0, 5.4d0, &
         &     2.0d0, 1.8d0, 1.0d0   ]

    ! U-238
    real(real64),   dimension(mq), parameter :: u238 = &
         & [   1.0d0, 1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0, &
         &     1.0d0, 0.989d0, 0.987d0, 0.985d0, 0.962d0, 0.960d0, 0.969d0, &
         &     0.970d0, 0.975d0, 0.978d0   ]
    real(real64),   dimension(mq), parameter :: u238f1 = &
         & [ 1.63d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0,  1.0d0, 1.0d0, &
         &   1.1d0,  1.4d0, 2.1d0, 3.6d0, 8.8d0, 10.4d0, 6.2d0, &
         &   2.8d0, 1.8d0, 1.0d0   ]

    ! Np-237
    real(real64),   dimension(mq), parameter :: np237 = &
         & [   1.0d0,   1.0d0,   1.0d0,   1.0d0,   0.998d0, 0.992d0, 0.990d0, &
         &     0.990d0, 0.940d0, 0.925d0, 0.938d0, 0.900d0, 0.900d0, 0.911d0, &
         &     0.923d0, 0.946d0, 0.966d0   ]
    real(real64),   dimension(mq), parameter :: np237f1 = &
         & [   3.4d0, 1.0d0, 1.0d0, 1.0d0,  1.2d0,  1.8d0, 1.8d0, &
         &     1.7d0, 2.1d0, 2.6d0, 6.3d0, 11.3d0, 11.2d0, 5.8d0, &
         &     1.9d0, 1.d0,  1.d0   ]

    ! Pu-239
    real(real64),   dimension(mq), parameter :: pu239 = &
         & [   1.0d0,   1.0d0,   1.0d0,  1.0d0,   0.960d0, 0.925d0, 0.910d0, &
         &     0.900d0, 0.955d0, 0.92d0, 0.915d0, 0.900d0, 0.900d0, 0.905d0, &
         &     0.915d0, 0.935d0, 0.930d0   ]
    real(real64),   dimension(mq), parameter :: pu239f1 = &
         & [   3.0d0, 1.0d0, 1.0d0,  1.0d0,  1.5d0,  2.5d0,  2.6d0, &
         &     2.6d0, 3.5d0, 6.2d0, 12.0d0, 25.0d0, 26.9d0, 17.3d0, &
         &     9.5d0, 7.0d0, 2.5d0   ]

! ======================================================================

    z0 = ztar0
    t0 = ener0
    if (z0 <= zfit(1))     then
       iz = 1
    elseif (z0 >= zfit(5)) then
       iz = 5
    else
       do i = 1,5
          if (z0 <= zfit(i)) then
             if ((z0-zfit(i-1)) < (zfit(i)-z0))  then
                iz = i - 1
             else
                iz = i
             endif
             go to 10
          endif
       enddo
    endif
10  continue
    ia = nint(atar0) 
    do  i = 1,mq
       xq(i) = tfit(i)
       if (iz == 1)     then
          yq(i) = ac227(i)
          zq(i) = ac227f1(i)
       elseif (iz == 2) then
          yq(i) = th232(i)
          zq(i) = th232f1(i)
       elseif (iz == 3) then
          if (ia <= 233)     then 
             yq(i) = u233(i) 
             zq(i) = u233f1(i) 
          elseif (ia == 235) then
             yq(i) = u235(i)
             zq(i) = u235f1(i)
          else
             yq(i) = u238(i)
             zq(i) = u238f1(i)
          endif
       elseif (iz == 4) then
          yq(i) = np237(i)
          zq(i) = np237f1(i)
       elseif (iz == 5) then
          yq(i) = pu239(i)
          zq(i) = pu239f1(i)
       endif
    enddo
    x = t0       
!   Exclude extrapolation:     
    if (x < xq(1)) then
       y = yq(1)
       z = zq(1)
       go to 30
    elseif (x > xq(mq)) then 
       y = yq(mq)
       z = zq(mq)
       go to 30
    endif

    do i = 1,mq 
       k = i 
       if (abs(x-xq(i)) < 1.d-10) then
          y = yq(k)
          z = zq(k) 
          go to 30
       endif
    end do
    k = 1 
20  if (x < xq(k)) then
       if (k <= 1) then
          y1 = yq(1) 
          z1 = zq(1)
          x1 = xq(1) 
          y2 = yq(2) 
          z2 = zq(2) 
          x2 = xq(2) 
          y3 = yq(3) 
          z3 = zq(3) 
          x3 = xq(3) 
       else
!         if (k-(mq-1))   15,14,14 
          if (k >= (mq-1)) then
             y1 = yq(mq-2) 
             z1 = zq(mq-2)
             x1 = xq(mq-2) 
             y2 = yq(mq-1) 
             z2 = zq(mq-1) 
             x2 = xq(mq-1) 
             y3 = yq(mq) 
             z3 = zq(mq) 
             x3 = xq(mq) 
          else
             y1 = yq(k-1) 
             z1 = zq(k-1)
             x1 = xq(k-1) 
             y2 = yq(k) 
             z2 = zq(k) 
             x2 = xq(k) 
             y3 = yq(k+1) 
             z3 = zq(k+1) 
             x3 = xq(k+1) 
          endif
       endif
    elseif (x == xq(k)) then
       y = yq(k)
       z = zq(k) 
       go to 30
    elseif (x > xq(k)) then
       k = k + 1 
       if (k > mq) then
          y1 = yq(mq-2) 
          z1 = zq(mq-2)
          x1 = xq(mq-2) 
          y2 = yq(mq-1) 
          z2 = zq(mq-1) 
          x2 = xq(mq-1) 
          y3 = yq(mq) 
          z3 = zq(mq) 
          x3 = xq(mq) 
       else
          go to 20 
       endif
    endif

    x23 = x2 - x3
    x31 = x3 - x1
    x12 = x1 - x2
    d  = x23*x1**2 + x31*x2**2 + x12*x3**2 
    da = x23*y1 + x31*y2 + x12*y3 
    db = (y2 - y3)*x1**2 + (y3 - y1)*x2**2 + (y1 - y2)*x3**2 
    dc = (x2*y3 - x3*y2)*x1**2 + (x3*y1 - x1*y3)*x2**2 + &
         &     (x1*y2 - x2*y1)*x3**2 
    da1 = x23*z1 + x31*z2 + x12*z3 
    db1 = (z2 - z3)*x1**2 + (z3 - z1)*x2**2 + (z1 - z2)*x3**2 
    dc1 = (x2*z3 - x3*z2)*x1**2 + (x3*z1 - x1*z3)*x2**2 + &
         &      (x1*z2 - x2*z1)*x3**2
    if (abs(d) < 1.0d-20) then
       write ( *, 2000) "258-264"
       d= 1.0d-20
    end if
    a = da/d
    b = db/d
    c = dc/d
    y = a*x**2 + b*x + c 
    a1 = da1/d
    b1 = db1/d
    c1 = dc1/d
    z = a1*x**2 + b1*x + c1 

30  continue
    fitaf = y 
    fitaf1= z 
    return 

! ======================================================================
2000 format(3x, "Warning: Divide by zero error prevented in ", &
          & "'fita.f90', line ", A, ". Fission interpolation error occurred.")
! ======================================================================
  end subroutine fitafac


  subroutine fitafpa (atar0, ztar0, ener0, fitaf, fitaf1)

! ======================================================================
!
!    Quadratic interpolation of fitaf = f(z0, t0) for subactinides.
!
!     Written by S. G. Mashnik
!    "Last" change: 12-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
!    Changed by K. K. Gudima, using new fit by M. I. Baznat, Nov., 2004
!    Edited by A. J. Sierk, LANL T-16, January, 2005.
!    Entries for W (182,184,186) added by Baznat and Gudima, June, 2005
!      (W183 modified.)
!    Error in entries for 181 Ta corrected by AJS, August 2011.
!    Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use evaporationParams, only: one

    implicit none
    real(real64), intent(in)  :: atar0
    real(real64), intent(in)  :: ztar0
    real(real64), intent(in)  :: ener0
    real(real64), intent(out) :: fitaf
    real(real64), intent(out) :: fitaf1

    integer(int32) :: ia, i, iz, k
    real(real64)   :: a, b, c, df21, df32, t0, x, x1, x2, x21, x3, x31, &
         & x32, y, y1, y2, y21, y3, y32, z0

! ======================================================================

    integer(int32), parameter :: mq = 21

    ! For calculation
    real(real64), dimension(mq) :: xq   = 0.d0
    real(real64), dimension(mq) :: yq   = 0.d0

    ! For fits
    real(real64), dimension(mq), parameter :: tfit = &
         & [   0.035d0, 0.04d0, 0.045d0, 0.05d0, 0.055d0, 0.06d0, 0.065d0, &
         &     0.07d0,  0.1d0,  0.15d0,  0.2d0,  0.3d0,   0.4d0,  0.5d0, &
         &     0.8d0,   1.0d0,  1.5d0,   2.0d0,  2.5d0,  3.0d0,  5.d0   ]
    real(real64), dimension(13), parameter :: zfit = &
         & [   67.d0, 70.d0, 73.d0, 74.d0, 75.d0, 78.d0, 79.d0, &
         &     80.d0, 81.d0, 82.d0, 83.d0, 84.d0, 85.d0   ]




    ! Ho-165
    real(real64), dimension(mq), parameter :: ho165 = &
         & [   1.004d0, 1.004d0, 1.004d0, 1.004d0, 1.004d0, 1.004d0, 1.004d0, &
         &     1.004d0, 1.004d0, 1.004d0, 0.985d0, 0.959d0, 0.949d0, 0.945d0, &
         &     0.948d0, 0.956d0, 0.989d0, 1.060d0, 1.163d0, 1.3d0,   1.3d0   ]
    ! Yb-173
    real(real64), dimension(mq), parameter :: yb173 = &
         & [   1.04d0,  1.04d0,  1.04d0,  1.04d0,  1.04d0,  1.04d0,  1.04d0, &
         &     1.04d0,  1.04d0,  1.04d0,  1.026d0, 0.993d0, 0.968d0, 0.953d0, &
         &     0.923d0, 0.917d0, 0.921d0, 0.932d0, 0.942d0, 0.953d0, 1.007d0   ]

    ! Ta-181
    real(real64), dimension(mq), parameter :: ta181 = &
!   Normal cem03.xx values; wrong!  AJS (6/20/11)
!         & [   1.059d0, 1.059d0, 1.059d0, 1.059d0, 1.059d0, 1.059d0, 1.059d0, &
!         &     1.059d0, 1.059d0, 1.059d0, 1.038d0, 0.992d0, 0.967d0, 0.943d0, &
!         &     0.897d0, 0.886d0, 0.877d0, 0.879d0, 0.886d0, 0.893d0, 0.919d0   ]
! (8/09/11) Modified to get Prokofiev systematics for 181 Ta (AJS):
         & [   1.000d0, 1.000d0,  1.000d0,  1.000d0,  1.005d0,  1.010d0,  1.015d0, &
         &     1.020d0, 1.071d0,  1.050d0,  1.0396d0, 1.0159d0, 0.9992d0, 0.9844d0, &
         &     0.949d0, 0.9349d0, 0.9221d0, 0.9219d0, 0.9248d0, 0.9288d0, 0.9432d0   ]

    ! W-182
    real(real64), dimension(mq), parameter :: w182 = &
         &  [   1.115d0, 1.108d0, 1.101d0, 1.096d0, 1.093d0, 1.090d0, 1.088d0, &
         &      1.086d0, 1.079d0, 1.053d0, 1.041d0, 1.020d0, 1.006d0, 0.995d0, &
         &      0.963d0, 0.949d0, 0.932d0, 0.929d0, 0.929d0, 0.930d0, 0.942d0   ]
    ! W-183
    real(real64), dimension(mq), parameter :: w183 = &
         &  [   1.095d0, 1.094d0, 1.091d0, 1.090d0, 1.089d0, 1.088d0, &
         &      1.087d0, 1.086d0, 1.069d0, 1.053d0, 1.037d0, 1.013d0, &
         &      1.001d0, 0.987d0, 0.955d0, 0.939d0, 0.923d0, 0.920d0, &
         &      0.921d0, 0.924d0, 0.936d0   ]
    ! W-184
    real(real64), dimension(mq), parameter :: w184 = &
         &  [   1.110d0, 1.105d0, 1.100d0, 1.090d0, 1.085d0, 1.082d0, &
         &      1.079d0, 1.076d0, 1.070d0, 1.053d0, 1.037d0, 1.014d0, &
         &      0.999d0, 0.9865d0,0.951d0, 0.935d0, 0.917d0, 0.9135d0, &
         &      0.914d0, 0.9165d0,0.933d0   ]
    ! W-186
    real(real64), dimension(mq), parameter :: w186 = &
         &  [   1.110d0, 1.105d0, 1.100d0, 1.090d0, 1.090d0, 1.084d0, &
         &      1.078d0, 1.075d0, 1.068d0, 1.047d0, 1.031d0, 1.006d0, &
         &      0.991d0, 0.976d0, 0.937d0, 0.919d0, 0.900d0, 0.900d0, &
         &      0.900d0, 0.902d0, 0.919d0   ]

    ! Re-186
    real(real64), dimension(mq), parameter :: re186 = &
         & [   1.027d0, 1.027d0, 1.027d0, 1.027d0, 1.027d0, 1.027d0, &
         &     1.027d0, 1.027d0, 1.027d0, 1.027d0, 1.025d0, 1.0015d0, &
         &     1.001d0, 0.993d0, 0.965d0, 0.950d0, 0.931d0, 0.915d0, &
         &     0.924d0, 0.924d0, 0.9355d0   ]

    ! Pt-195
    real(real64), dimension(mq), parameter :: pt195 = &
         & [   1.021d0, 1.021d0, 1.021d0,  1.021d0, 1.021d0, 1.021d0, &
         &     1.021d0, 1.021d0, 1.021d0,  1.021d0, 1.020d0, 1.009d0, &
         &     1.002d0, 0.997d0, 0.9794d0, 0.970d0, 0.950d0, 0.939d0, &
         &     0.935d0, 0.934d0, 0.939d0   ]

    ! Au-197
    real(real64), dimension(mq), parameter :: au197 = &
         & [   1.105d0,  1.086d0,  1.055d0, 1.067d0, 1.065d0, 1.066d0, &
         &     1.061d0,  1.054d0,  1.036d0, 1.024d0, 1.022d0, 1.0105d0, &
         &     1.006d0,  0.9995d0, 0.983d0, 0.975d0, 0.958d0, 0.948d0, &
         &     0.9443d0, 0.943d0,  0.946d0   ]

    ! Hg-202
    real(real64), dimension(mq), parameter :: hg202 = &
         & [   1.075d0, 1.062d0,  1.059d0, 1.062d0, 1.069d0, 1.068d0, &
         &     1.066d0, 1.059d0,  1.045d0, 1.027d0, 1.019d0, 1.006d0, &
         &     1.001d0, 0.997d0,  0.982d0, 0.974d0, 0.958d0, 0.948d0, &
         &     0.942d0, 0.9404d0, 0.943d0   ]

    ! Tl-205
    real(real64), dimension(mq), parameter :: tl205 = &
         & [   1.005d0, 1.031d0, 1.032d0,  1.0335d0, 1.043d0, 1.044d0, &
         &     1.041d0, 1.035d0, 1.019d0,  1.001d0,  0.991d0, 0.983d0, &
         &     0.982d0, 0.978d0, 0.9705d0, 0.966d0,  0.955d0, 0.946d0, &
         &     0.940d0, 0.938d0, 0.9405d0   ]

    ! Pb-204
    real(real64), dimension(mq), parameter :: pb204 = &
         & [   1.018d0,  1.015d0,  1.011d0, 1.014d0, 1.017d0, 1.020d0, &
         &     1.017d0,  1.014d0,  1.007d0, 1.007d0, 0.999d0, 0.989d0, &
         &     0.9867d0, 0.982d0,  0.976d0, 0.975d0, 0.969d0, 0.964d0, &
         &     0.960d0,  0.9604d0, 0.969d0   ]
    ! Pb-206
    real(real64), dimension(mq), parameter :: pb206 = &
         & [   1.015d0,  1.017d0, 1.017d0,  1.021d0,  1.026d0, 1.025d0, &
         &     1.022d0,  1.018d0, 1.008d0,  1.000d0,  0.993d0, 0.981d0, &
         &     0.979d0,  0.975d0, 0.9685d0, 0.9685d0, 0.962d0, 0.956d0, &
         &     0.9525d0, 0.951d0, 0.955d0   ]
    ! Pb-207
    real(real64), dimension(mq), parameter :: pb207 = &
         & [   1.008d0,  1.007d0,  1.016d0, 1.021d0, 1.029d0, 1.029d0, &
         &     1.024d0,  1.019d0,  1.009d0, 0.998d0, 0.988d0, 0.9755d0, &
         &     0.972d0,  0.968d0,  0.962d0, 0.962d0, 0.956d0, 0.950d0, &
         &     0.9446d0, 0.9435d0, 0.945d0   ]
    ! pb-208
    real(real64), dimension(mq), parameter :: pb208 = &
         & [   1.005d0, 1.022d0, 1.028d0, 1.025d0, 1.033d0, 1.033d0, &
         &     1.028d0, 1.028d0, 1.013d0, 1.000d0, 0.992d0, 0.978d0, &
         &     0.974d0, 0.970d0, 0.964d0, 0.962d0, 0.956d0, 0.949d0, &
         &     0.945d0, 0.942d0, 0.945d0   ]

    ! Bi-209
    real(real64), dimension(mq), parameter :: bi209 = &
         & [   0.999d0, 1.004d0, 1.006d0, 1.010d0, 1.016d0, 1.017d0, &
         &     1.013d0, 1.010d0, 0.996d0, 0.987d0, 0.979d0, 0.968d0, &
         &     0.963d0, 0.958d0, 0.956d0, 0.957d0, 0.959d0, 0.955d0,  &
         &     0.952d0, 0.949d0, 0.953d0   ]

    ! Po-210
    real(real64), dimension(mq), parameter :: po210 = &
         & [   1.039d0, 1.039d0, 1.039d0,  1.039d0,  1.039d0,  1.039d0, &
         &     1.039d0, 1.046d0, 1.043d0,  1.020d0,  1.001d0,  0.976d0, &
         &     0.969d0, 0.965d0, 0.9665d0, 0.9695d0, 0.9743d0, 0.9727d0, &
         &     0.973d0, 0.973d0, 0.992d0   ]

    ! At-211
    real(real64), dimension(mq), parameter :: at211 = &
         & [   1.700d0, 1.700d0, 1.700d0, 1.700d0, 1.700d0, 1.700d0, &
         &     1.700d0, 1.700d0, 1.700d0, 1.550d0, 1.445d0, 1.312d0, &
         &     1.243d0, 1.188d0, 1.107d0, 1.083d0, 1.055d0, 1.050d0, &
         &     1.060d0, 1.075d0, 1.205d0   ]

! ======================================================================

    z0 = ztar0
    t0 = ener0
    if (z0 <= zfit(1))      then
       iz = 1
    elseif (z0 >= zfit(13))  then
       iz = 13
    else
       do i = 1,13
          if (z0 <= zfit(i)) then
             if ((z0-zfit(i-1)) < (zfit(i)-z0))  then
                iz = i - 1
             else
                iz = i
             endif
             go to 10
          endif
       enddo
    endif
10  continue
    ia = nint(atar0)  
    do  i = 1,mq
       xq(i) = tfit(i)
       if (iz == 1)     then
          yq(i) = ho165(i)
       elseif (iz == 2) then
          yq(i) = yb173(i)
       elseif (iz == 3) then
          yq(i) = ta181(i)
       elseif (iz == 4) then
          if (ia <= 182) then
             yq(i) =  w182(i)
          elseif (ia == 183) then
             yq(i) =  w183(i)
          elseif (ia == 184) then
             yq(i) =  w184(i)
          else
             yq(i) =  w186(i)
          endif
       elseif (iz == 5) then
          yq(i) = re186(i)
       elseif (iz == 6) then
          yq(i) = pt195(i)
       elseif (iz == 7) then
          yq(i) = au197(i)
       elseif (iz == 8) then
          yq(i) = hg202(i)
       elseif (iz == 9) then
          yq(i) = tl205(i)
       elseif (iz == 10) then
          if (ia <= 204) then 
             yq(i) = pb204(i)
          elseif (ia == 206) then 
             yq(i) = pb206(i)
          elseif (ia == 207) then 
             yq(i) = pb207(i)
          else
             yq(i) = pb208(i)
          endif
       elseif (iz == 11) then
          yq(i) = bi209(i)
       elseif (iz == 12) then
          yq(i) = po210(i)
       elseif (iz == 13) then
          yq(i) = at211(i)
       endif
    enddo
    x = t0  
!   Exclude extrapolation:     
    if (x < xq(1))  then
       y = yq(1)
       go to 30
    endif
    if (x > xq(mq)) then 
       y = yq(mq)
       go to 30
    endif

    do i = 1,mq 
       k = i 
       if (abs(x - xq(i)) < 1.d-10) then
          y = yq(k)
          go to 30
       endif
    end do

    k = 1 
20  if (x < xq(k)) then
       if (k <= 1) then
          y1 = yq(1) 
          x1 = xq(1) 
          y2 = yq(2) 
          x2 = xq(2) 
          y3 = yq(3) 
          x3 = xq(3) 
       else
          if (k >= (mq - 1)) then
             y1 = yq(mq-2) 
             x1 = xq(mq-2) 
             y2 = yq(mq-1) 
             x2 = xq(mq-1) 
             y3 = yq(mq) 
             x3 = xq(mq) 
          else
             y1 = yq(k-1) 
             x1 = xq(k-1) 
             y2 = yq(k) 
             x2 = xq(k) 
             y3 = yq(k+1) 
             x3 = xq(k+1) 
          endif
       endif
    elseif (x == xq(k)) then
       y = yq(k) 
       go to 30
    else
       k = k + 1 
       if (k > mq) then
          y1 = yq(mq-2) 
          x1 = xq(mq-2) 
          y2 = yq(mq-1) 
          x2 = xq(mq-1) 
          y3 = yq(mq) 
          x3 = xq(mq) 
       else
          go to 20 
       endif
    endif
!  Modified 7/19/11 by AJS:
    x32 = x3 - x2
    x31 = x3 - x1
    x21 = x2 - x1
    y32 =  y3 -  y2
    y21 =  y2 -  y1
    if (x21 < 1.0d-20 .and. x21 > -1.0d-20) then
       x21 = 1.0d-20
       write(*,2000) '596'
    end if
    if (x31 < 1.0d-20 .and. x31 > -1.0d-20) then
       x31 = 1.0d-20
       write(*,2000) '600'
    end if
    if (x32 < 1.0d-20 .and. x32 > -1.0d-20) then
       x32 = 1.0d-20
       write(*,2000) '604'
    end if
    df32 = y32/x32 
    df21 = y21/x21 
    a = (df32 - df21)/x31
    b = (df21*(x2 + x3) - df32*(x1 + x2))/x31
    c = y1*x2*x3/(x31*x21) - y2*x1*x3/(x21*x32) + &
         &    y3*x1*x2/(x31*x32)

    y = a*x**2 + b*x + c

30  continue 
    fitaf = y 
    fitaf1 = one
    return 

! ======================================================================
2000 format(3x, "Warning: Divide by zero error prevented in ", &
          & "'fita.f90', line ", A, ". Fission interpolation error occurred.")
! ======================================================================
  end subroutine fitafpa
      
         
! ======================================================================

  SUBROUTINE FITAFPAQ(atar0,ztar0,ener0,fitaf,fitaf1)

! ======================================================================
!
!      quadratic interpolation of fitaf=f(Z0,T0) for subactinides
!         (used with the Modified DCM/QGSM hybrid)
!
!
! Written by KKG, 12/12/05
! Modified by CMJ, XCP-3, July 2018 (Evaporation class creation)
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64

    implicit none
    real(real64), intent(in)  :: atar0
    real(real64), intent(in)  :: ztar0
    real(real64), intent(in)  :: ener0
    real(real64), intent(out) :: fitaf
    real(real64), intent(out) :: fitaf1

    integer(int32) :: i, ia, iz, k
    real(real64)   :: a, b, c, d, da, db, dc, t0, x, x1, x2, x3, y, &
         & y1, y2, y3, z0

! ======================================================================

    integer(int32), parameter :: mq = 22

    ! For calculation
    real(real64), dimension(mq) :: xq = 0.d0
    real(real64), dimension(mq) :: yq = 0.d0

    ! For fits
    real(real64), dimension(mq), parameter :: tfit = &
         & [    0.035, 0.040, 0.045, 0.050, 0.055, 0.060, &
         &      0.065, 0.070, 0.100, 0.150, 0.200, 0.300, &
         &      0.400, 0.500, 0.800, 1.000, 1.500, 2.000, &
         &      2.500, 3.000, 5.000, 10.000   ]
    real(real64), dimension(13), parameter :: zfit = &
         & [   67.0, 70.0, 73.0, 74.0, 75.0, 78.0, &
         &     79.0, 80.0, 81.0, 82.0, 83.0, 84.0, &
         &     85.0   ]


    ! Ho-165
    real(real64), dimension(mq), parameter :: ho165 = &
         & [   1.040, 1.040, 1.040, 1.040, 1.040, 1.040, &
         &     1.040, 1.040, 1.040, 1.040, 1.009, 0.987, &
         &     0.986, 0.983, 0.994, 1.004, 1.055, 1.135, &
         &     1.250, 1.306, 1.700, 2.500   ]

    ! Yb-173
    real(real64), dimension(mq), parameter :: yb173 = &
         & [   1.070, 1.070,  1.070, 1.070, 1.070, 1.070, &
         &     1.070, 1.070,  1.070, 1.070, 1.045, 1.015, &
         &     0.994, 0.9735, 0.945, 0.940, 0.941, 0.956, &
         &     0.968, 0.979,  1.021, 1.065   ]

    ! Ta-181
    real(real64), dimension(mq), parameter :: ta181 = &
         & [   1.086, 1.086, 1.086, 1.086, 1.086,  1.086, &
         &     1.086, 1.086, 1.086, 1.086, 1.0645, 1.04,  &
         &     1.023, 1.01,  0.969, 0.950, 0.931,  0.927, &
         &     0.926, 0.929, 0.945, 0.960   ]

    ! W-183
    real(real64), dimension(mq), parameter :: w183 = &
         & [   1.069, 1.069, 1.069, 1.069, 1.069, 1.069,  &
         &     1.069, 1.069, 1.069, 1.069, 1.053, 1.0335, &
         &     1.021, 1.009, 0.935, 0.959, 0.935, 0.926,  &
         &     0.925, 0.925, 0.940, 0.952   ]

    ! Re-186
    real(real64), dimension(mq), parameter :: re186 = &
         & [   1.050, 1.050, 1.050, 1.050, 1.050, 1.050, &
         &     1.050, 1.050, 1.050, 1.050, 1.044, 1.031, &
         &     1.024, 1.015, 0.987, 0.974, 0.946, 0.931, &
         &     0.926, 0.923, 0.934, 0.950   ]

    ! Pt-195
    real(real64), dimension(mq), parameter :: pt195 = &
         & [   1.038,  1.038, 1.038, 1.038, 1.038, 1.038, &
         &     1.038,  1.038, 1.038, 1.038, 1.040, 1.032, &
         &     1.0275, 1.022, 1.004, 0.994, 0.972, 0.955, &
         &     0.943,  0.938, 0.938, 0.953   ]

    ! Au-197
    real(real64), dimension(mq), parameter :: au197 = &
         & [   1.142, 1.113,  1.102,  1.108, 1.103, 1.095, &
         &     1.091, 1.082,  1.060,  1.044, 1.044, 1.036, &
         &     1.031, 1.0237, 1.007,  0.999, 0.981, 0.967, &
         &     0.957, 0.951,  0.9515, 0.967   ]

    ! Hg-202
    real(real64), dimension(mq), parameter :: hg202 = &
         & [   1.087, 1.0955, 1.093, 1.095, 1.098, 1.097, &
         &     1.093, 1.084,  1.06,  1.045, 1.044, 1.033, &
         &     1.029, 1.023,  1.005, 0.998, 0.982, 0.968, &
         &     0.956, 0.948,  0.945, 0.957   ]

    ! Tl-205
    real(real64), dimension(mq), parameter :: tl205 = &
         & [   1.077,  1.055, 1.064, 1.076, 1.072, 1.072, &
         &     1.067,  1.061, 1.038, 1.02,  1.015, 1.012, &
         &     1.0085, 1.005, 0.994, 0.988, 0.978, 0.966, &
         &     0.956,  0.948, 0.944, 0.953   ]

    ! Pb-204
    real(real64), dimension(mq), parameter :: pb204 = &
         & [   1.022, 1.022, 1.024, 1.028, 1.03,   1.032, &
         &     1.03,  1.027, 1.018, 1.023, 1.024,  1.016, &
         &     1.011, 1.004, 0.997, 0.994, 0.9885, 0.984, &
         &     0.977, 0.974, 0.973, 0.984   ]
    ! Pb-206
    real(real64), dimension(mq), parameter :: pb206 = &
         & [   1.02,  1.026, 1.03,  1.029, 1.038, 1.036, &
         &     1.034, 1.032, 1.020, 1.014, 1.016, 1.009, &
         &     1.004, 0.998, 0.991, 0.988, 0.983, 0.976, &
         &     0.970, 0.965, 0.961, 0.971   ]
    ! Pb-207
    real(real64), dimension(mq), parameter :: pb207 = &
         & [   1.005, 1.0275, 1.03,  1.037,  1.044, 1.045, &
         &     1.042, 1.038,  1.022, 1.0155, 1.01,  1.002, &
         &     0.998, 0.993,  0.985, 0.982,  0.977, 0.970, &
         &     0.962, 0.957,  0.951, 0.960   ]
    ! Pb-208
    real(real64), dimension(mq), parameter :: pb208 = &
         & [   1.036, 1.028, 1.032, 1.039, 1.043,  1.048, &
         &     1.045, 1.04,  1.027, 1.014, 1.012,  1.006, &
         &     1.001, 0.995, 0.988, 0.984, 0.9778, 0.969, &
         &     0.961, 0.955, 0.949, 0.958   ]

    ! Bi-209
    real(real64), dimension(mq), parameter :: bi209 = &
         & [   1.017, 1.018,  1.029,  1.03,  1.033, 1.033, &
         &     1.031, 1.028,  1.011,  1.006, 1.009, 0.999, &
         &     0.991, 0.984,  0.9775, 0.978, 0.977, 0.973, &
         &     0.968, 0.9635, 0.960,  0.967   ]

    ! Po-210
    real(real64), dimension(mq), parameter :: po210 = &
         & [    1.055, 1.055, 1.055, 1.055, 1.055, 1.055, &
         &      1.055, 1.055, 1.06,  1.045, 1.035, 1.011, &
         &      0.996, 0.991, 0.985, 0.986, 0.988, 0.989, &
         &      0.986, 0.983, 0.984, 0.996   ]

    ! At-211
    real(real64), dimension(mq), parameter :: at211 = &
         & [   1.27,  1.27,   1.27,  1.27,  1.27,  1.27,  &
         &     1.27,  1.27,   1.220, 1.152, 1.108, 1.062, &
         &     1.033, 1.018,  1.01,  1.011, 1.015, 1.016, &
         &     1.0145, 1.013, 1.018, 1.045   ]

! ======================================================================

    Z0=ztar0
    T0=ener0
    if(Z0 <= ZFIT(1))      then
       IZ=1
    elseif(Z0 >= ZFIT(13))  then
       IZ=13
    else
       do  i=1,13
          if(Z0 <= ZFIT(i)) then
             if((Z0-ZFIT(i-1)) < (ZFIT(i)-Z0))  then
                IZ=i-1
             else
                IZ=i
             endif
             go  to  1
          endif
       enddo
    endif
1   continue   
    IA=INT(atar0+0.1)        
    do  i=1,MQ
       XQ(i)=TFIT(i)
       if(IZ == 1)        then
          YQ(i)=Ho165(i)
       elseif(IZ == 2)    then
          YQ(i)=Yb173(i)
       elseif(IZ == 3)    then
          YQ(i)=Ta181(i)
       elseif(IZ == 4)    then
          YQ(i)= W183(i)
       elseif(IZ == 5)    then
          YQ(i)=Re186(i)
       elseif(IZ == 6)    then
          YQ(i)=Pt195(i)
       elseif(IZ == 7)    then
          YQ(i)=Au197(i)
       elseif(IZ == 8)    then
          YQ(i)=Hg202(i)
       elseif(IZ == 9)    then
          YQ(i)=Tl205(i)
       elseif(IZ == 10)   then
          if(IA <= 204)      then 
             YQ(i)=Pb204(i)
          elseif(IA == 206)  then 
             YQ(i)=Pb206(i)
          elseif(IA == 207)  then 
             YQ(i)=Pb207(i)
          else
             YQ(i)=Pb208(i)
          endif
       elseif(IZ == 11)    then
          YQ(i)=Bi209(i)
       elseif(IZ == 12)    then
          YQ(i)=Po210(i)
       elseif(IZ == 13)    then
          YQ(i)=At211(i)
       endif
    enddo
    X=T0  
!      exclude extrapolation!                       
    if(X < XQ(1))     then
       Y=YQ(1)
       go  to  19 
    endif
    if(X.GT.XQ(MQ))    then 
       Y=YQ(MQ)
       go  to  19 
    endif
    do  i=1,MQ 
       K=i 
       if(ABS(X-XQ(i)) < 1.D-10) go  to  16 
    end do
    K=1 
10  IF( X-XQ(K) < 0 ) then
       go to 11
    else if ( X-XQ(K) == 0 ) then
       go to 16
    else 
       go to 17
    end if
11  IF(K-1 <= 0 ) then
       go to 12
    else
       go to 13
    end if
12  Y1=YQ(1) 
    X1=XQ(1) 
    Y2=YQ(2) 
    X2=XQ(2) 
    Y3=YQ(3) 
    X3=XQ(3) 
    GO  TO  18 
13  IF( K-(MQ-1) < 0 ) then
       go to 15
    else
       go to 14 
    end if
14  Y1=YQ(MQ-2) 
    X1=XQ(MQ-2) 
    Y2=YQ(MQ-1) 
    X2=XQ(MQ-1) 
    Y3=YQ(MQ) 
    X3=XQ(MQ) 
    GO  TO  18 
15  Y1=YQ(K-1) 
    X1=XQ(K-1) 
    Y2=YQ(K) 
    X2=XQ(K) 
    Y3=YQ(K+1) 
    X3=XQ(K+1) 
    GO  TO  18 
16  Y = YQ(K) 
    GO  TO  19 
17  K = K + 1 
    IF(K.GT.MQ)  GO  TO  14 
    GO  TO  10 
18  D =(X2-X3)*X1**2+(X3-X1)*X2**2+(X1-X2)*X3**2 
    DA=(X2-X3)*Y1   +(X3-X1)*Y2   +(X1-X2)*Y3 
    DB=(Y2-Y3)*X1**2+(Y3-Y1)*X2**2+(Y1-Y2)*X3**2 
    DC=(X2*Y3-X3*Y2)*X1**2+(X3*Y1-X1*Y3)*X2**2+(X1*Y2-X2*Y1)*X3**2 
    IF(ABS(D).LT.1.0D-20) then
       write(*,2000) '887'
       D = 1.0d-20
    end if
    A=DA/D 
    B=DB/D 
    C=DC/D 
    Y = A*X**2 + B*X + C
19  CONTINUE 
    fitaf=Y 
    fitaf1=1.
    RETURN 
! ======================================================================
2000 format(3x, "Warning: Divide by zero error prevented in ", &
          & "'fita.f90', line ", A, ". Fission interpolation error occurred.")
! ======================================================================
  END SUBROUTINE FITAFPAQ


! ======================================================================

  SUBROUTINE FITAFACQ(atar0,ztar0,ener0,fitaf,fitaf1)

! ======================================================================
!
!      quadratic interpolation of fitaf=f(Z0,T0) for actinides
!         (used with the Modified DCM/QGSM hybrid)
!
!
! Written by KKG, 12/12/05
! Modified by CMJ, XCP-3, July 2018 (Evaporation class creation)
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64

    implicit none
    real(real64), intent(in)  :: atar0
    real(real64), intent(in)  :: ztar0
    real(real64), intent(in)  :: ener0
    real(real64), intent(out) :: fitaf
    real(real64), intent(out) :: fitaf1

    integer(int32) :: i, ia, iz, k
    real(real64)   :: a, a1, b, b1, c, c1, d, da, da1, db, db1, dc, &
         & dc1, t0, x, x1, x2, x3, y, y1, y2, y3, z, z0, z1, z2, z3

! ======================================================================

    integer(int32), parameter :: mq = 18

    ! For calculation
    real(real64), dimension(mq) :: xq = 0.d0
    real(real64), dimension(mq) :: yq = 0.d0
    real(real64), dimension(mq) :: zq= 0.d0

    ! For fits
    real(real64), dimension(mq), parameter :: tfit = &
         & [   0.03, 0.04, 0.05, 0.07, 0.1,  0.15, &
         &     0.2,  0.3,  0.4,  0.5,  0.8,  1.0, &
         &     1.5,  2.0,  2.5,  3.0,  5.0, 10.0   ]
    real(real64), dimension( 5), parameter :: zfit = &
    & [   89.0, 90.0, 92.0, 93.0, 94.0   ]


    ! Ac-227
    real(real64), dimension(mq), parameter :: ac227 = &
         & [   1.150, 1.150, 1.150, 1.150, 1.150, 1.150, &
         &     1.12,  1.10,  1.025, 0.998, 0.983, 0.980, &
         &     0.983, 0.994, 0.997, 1.000, 1.003, 1.020   ]
    ! Ac-227f1
    real(real64), dimension(mq), parameter :: ac227f1 = &
         & [   1.00, 1.00, 1.00, 1.00, 1.00, 1.00, &
         &     1.00, 1.00, 1.00, 1.20, 5.70, 6.00, &
         &     1.00, 1.00, 1.00, 1.00, 1.00, 1.00   ]

    ! Th-232
    real(real64), dimension(mq), parameter :: th232 = &
         & [   1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  &
         &     1.00,  0.995, 0.991, 0.983, 0.976, 0.966, &
         &     0.955, 0.971, 0.974, 0.977, 0.979, 0.983   ]
    ! Th-232f1
    real(real64), dimension(mq), parameter :: th232f1 = &
         & [   1.00, 1.00, 1.00, 1.00,  1.00, 1.00, &
         &     1.00, 1.30, 2.00, 3.20, 13.6, 14.5, &
         &     9.30, 4.6,  1.9,  1.5,   1.4,  1.1   ]

    ! U233
    real(real64), dimension(mq), parameter :: u233 = &
         & [ 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, &
         &   1.000, 1.000, 0.990, 0.978, 0.973, 0.965, &
         &   0.96,  0.962, 0.977, 0.996, 1.008, 1.038   ]
    ! U233f1
    real(real64), dimension(mq), parameter :: u233f1 = &
         & [   1.00, 1.00, 1.00, 1.00,  1.00, 1.00, &
         &     1.00, 1.00, 2.00, 4.10, 10.1, 10.5,  &
         &     6.70, 2.00, 1.00, 1.00, 1.00, 1.00   ]

    ! U-235
    real(real64), dimension(mq), parameter :: u235 = &
         & [   1.000, 1.000, 1.000, 1.000, 1.000, 1.000, &
         &     1.000, 1.000, 0.995, 0.98,  0.975, 0.974, &
         &     0.964, 0.97,  0.983, 0.992, 1.010, 1.042   ]
   ! U-235f1
    real(real64), dimension(mq), parameter :: u235f1 = &
         & [   1.00, 1.00, 1.00, 1.00, 1.00, 1.00, &
         &     1.00, 1.00, 1.40, 2.90, 5.80, 7.60, &
         &     5.20, 1.90, 1.30, 1.00, 1.00, 1.00   ]

    ! U-238
    real(real64), dimension(mq), parameter :: u238 = &
         & [   1.000, 1.000, 1.000, 1.000, 1.000, 1.000, &
         &     1.000, 1.000, 1.000, 0.979, 0.975, 0.97,  &
         &     0.965, 0.965, 0.975, 0.985, 0.990, 1.010   ]
   ! U-238f1
    real(real64), dimension(mq), parameter :: u238f1 = &
         & [   1.00, 1.00, 1.00, 1.00, 1.00, 1.00, &
         &     1.00, 1.00, 1.00, 1.70, 4.40, 6.00, &
         &     4.70, 2.00, 1.60, 1.20, 1.00, 1.00   ]

    ! Np-237
    real(real64), dimension(mq), parameter :: np237 = &
         & [   1.000, 1.000, 1.000, 1.000, 1.000, 1.000, &
         &     1.000, 1.000, 0.99,  0.98,  0.963, 0.953, &
         &     0.952, 0.963, 0.965, 0.973, 0.991, 1.023   ]
    ! Np-237f1
    real(real64), dimension(mq), parameter :: np237f1 = &
         & [   1.00, 1.00, 1.00, 1.00,  1.00, 1.00, &
         &     1.00, 1.00, 1.70, 4.00, 10.5, 11.8,  &
         &    10.5,  5.8,  2.5,  1.00, 1.00, 1.00   ]

    ! Pu-239
    real(real64), dimension(mq), parameter :: pu239 = &
         & [   1.000, 1.000, 1.000, 1.000, 1.000, 1.000, &
         &     1.000, 1.000, 0.990, 0.962, 0.900, 0.900, &
         &     0.900, 0.918, 0.928, 0.946, 0.967, 0.988   ]
    ! Pu-239f1
    real(real64), dimension(mq), parameter :: pu239f1 = &
         & [   1.00, 1.00, 1.00, 1.00,  1.00, 1.00, &
         &     1.00, 1.00, 3.20, 7.70, 16.5, 19.0,  &
         &    14.0,  9.1,  4.40, 2.5,  1.00, 1.00   ]

! ======================================================================

    Z0=ztar0
    T0=ener0
    if(Z0 <= ZFIT(1))      then
       IZ=1
    elseif(Z0 >= ZFIT(5))  then
       IZ=5
    else
       do  i=1,5
          if(Z0 <= ZFIT(i)) then
             if((Z0-ZFIT(i-1)) < (ZFIT(i)-Z0))  then
                IZ=i-1
             else
                IZ=i
             endif
             go  to  1
          endif
       enddo
    endif
1   continue
    IA=INT(atar0+0.1)           
    do  i=1,MQ
       XQ(i)=TFIT(i)
       if(IZ == 1)        then
          YQ(i)=Ac227(i)
          ZQ(i)=Ac227f1(i)
       elseif(IZ == 2)    then
          YQ(i)=Th232(i)
          ZQ(i)=Th232f1(i)
       elseif(IZ == 3)    then
          if(IA <= 233)       then 
             YQ(i)=U233(i) 
             ZQ(i)=U233f1(i) 
          elseif(IA == 235)   then
             YQ(i)=U235(i)
             ZQ(i)=U235f1(i)
          else
             YQ(i)=U238(i)
             ZQ(i)=U238f1(i)
          endif
       elseif(IZ == 4)    then
          YQ(i)=Np237(i)
          ZQ(i)=Np237f1(i)
       elseif(IZ == 5)    then
          YQ(i)=Pu239(i)
          ZQ(i)=Pu239f1(i)
       endif
    enddo
    X=T0                      
!      exclude extrapolation!                       
    if(X < XQ(1))     then
       Y=YQ(1)
       Z=ZQ(1)
       go  to  19 
    endif
    if(X.GT.XQ(MQ))    then 
       Y=YQ(MQ)
       Z=ZQ(MQ)
       go  to  19 
    endif

    do  i=1,MQ 
       K=i 
       if(ABS(X-XQ(i)) < 1.D-10) go  to  16 
    end do
    K=1 
10  IF( X-XQ(K) < 0 ) then
       go to 11
    else if ( X-XQ(K) == 0 ) then
       go to 16
    else
       go to 17 
    end if
! -------------------------
11  IF(K-1 <= 0 ) then
       go to 12
    else
       go to 13
    end if
12  Y1=YQ(1) 
    Z1=ZQ(1)
    X1=XQ(1) 
    Y2=YQ(2) 
    Z2=ZQ(2) 
    X2=XQ(2) 
    Y3=YQ(3) 
    Z3=ZQ(3) 
    X3=XQ(3) 
    GO  TO  18 
13  IF(K-(MQ-1) < 0 ) then
       go to 15
    else
       go to 14 
    end if
14  Y1=YQ(MQ-2) 
    Z1=ZQ(MQ-2)
    X1=XQ(MQ-2) 
    Y2=YQ(MQ-1) 
    Z2=ZQ(MQ-1) 
    X2=XQ(MQ-1) 
    Y3=YQ(MQ) 
    Z3=ZQ(MQ) 
    X3=XQ(MQ) 
    GO  TO  18 
15  Y1=YQ(K-1) 
    Z1=ZQ(K-1)
    X1=XQ(K-1) 
    Y2=YQ(K) 
    Z2=ZQ(K) 
    X2=XQ(K) 
    Y3=YQ(K+1) 
    Z3=ZQ(K+1) 
    X3=XQ(K+1) 
    GO  TO  18 
16  Y = YQ(K)
    Z = ZQ(K) 
    GO  TO  19 
17  K = K + 1 
    IF(K.GT.MQ)  GO  TO  14 
    GO  TO  10 
18  D =(X2-X3)*X1**2+(X3-X1)*X2**2+(X1-X2)*X3**2 
    DA=(X2-X3)*Y1   +(X3-X1)*Y2   +(X1-X2)*Y3 
    DB=(Y2-Y3)*X1**2+(Y3-Y1)*X2**2+(Y1-Y2)*X3**2 
    DC=(X2*Y3-X3*Y2)*X1**2+(X3*Y1-X1*Y3)*X2**2+(X1*Y2-X2*Y1)*X3**2 
    DA1=(X2-X3)*Z1   +(X3-X1)*Z2   +(X1-X2)*Z3 
    DB1=(Z2-Z3)*X1**2+(Z3-Z1)*X2**2+(Z1-Z2)*X3**2 
    DC1=(X2*Z3-X3*Z2)*X1**2+(X3*Z1-X1*Z3)*X2**2+(X1*Z2-X2*Z1)*X3**2 
    IF(ABS(D).LT.1.0D-20) then
       write(*,2000) '1140'
       D = 1.0d-20
    end if
    A=DA/D 
    B=DB/D 
    C=DC/D 
    Y = A*X**2 + B*X + C 
    A1=DA1/D 
    B1=DB1/D 
    C1=DC1/D 
    Z = A1*X**2 + B1*X + C1 

19  CONTINUE 
    fitaf =Y 
    fitaf1=Z 
    if(fitaf1 < 1.0) fitaf1=1.0

    RETURN 
! ======================================================================
2000 format(3x, "Warning: Divide by zero error prevented in ", &
          & "'fita.f90', line ", A, ". Fission interpolation error occurred.")
! ======================================================================
  END SUBROUTINE FITAFACQ
