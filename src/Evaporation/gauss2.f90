
  function gauss2 (evapObj, xmean, sd)

! ======================================================================
!
!   This routine finds a random deviate with a gaussian distribution
!   of mean xmean and standard deviation sd.
!
!   This is a modified version of GASDEV, from Numerical Recipes, 2nd
!   Ed. by Press et al., Oxford (1992).  This version uses the function
!///RNDM as defined by Mokhov, instead of RAN1 from Numer. Recipes.////
!   RANG as used in MCNP6.
!
!   GAUSS2 uses an average of 4/pi = 1.273... random numbers for each
!   gaussian-distributed random number, as opposed to 12 for GAUSSN!!
!
!   Written by A. J. Sierk, LANL T-16, September, 2003.
!   Modified by A. J. S. December, 2011.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use evaporationParams, only: zro, one, two

    implicit none
    class(Evaporation), intent(inout) :: evapObj
    real(real64), intent(in   ) :: xmean   ! Mean value         of Gaussian
    real(real64), intent(in   ) :: sd      ! Standard Deviation of Gaussian
    real(real64)                :: gauss2  ! Random deviate of the Gaussian

    real(real64)   :: fac, gasdev, rsq, v1, v2

! ======================================================================

    integer(int32), save :: iset = 0_int32   ! Flags first or second call to the function
    real(real64),   save :: gset = zro       ! Contains random deviate for the second call to the function

! ======================================================================

    if (iset == 0) then
       ! Obtain valid 'rsq' value [in range (0, 1) ]
10     continue
       v1 = two*evapObj%rang() - one
       v2 = two*evapObj%rang() - one
       rsq = v1**2 + v2**2
       if (rsq >= one .or. rsq == zro) go to 10

       ! Obtain random deviate
       fac = sqrt(-two*log(rsq)/rsq)
       gset = v1*fac     ! Random deviate for second call
       gasdev = v2*fac   ! Random deviate for first call
       iset = 1          ! Flag to use "second" call
    else
       gasdev = gset   ! Set random deviate from previous call
       iset = 0        ! Reset to "first" call
    endif

    ! Obtain Value of gauss2 given the random deviate, mean, and standard deviation
    gauss2 = sd*gasdev + xmean

    return
! ======================================================================
  end function gauss2
