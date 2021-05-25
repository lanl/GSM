! copyright LANS/LANL/DOE - see file COPYRIGHT_INFO

module coalescenceParams
  !=======================================================================================
  ! Description:
  ! Various parameters used by many routines in the Coalescence class.
  ! ALL PARAMETERS IN THIS MODULE ARE PUBLIC.
  ! DO NOT PUT ANY FUNCTIONS OR SUBROUTINES IN THIS MODULE.
  !
  ! Organization:
  !   * Kind parameters
  !   * Numbers
  !   * General code & control
  !   * Parallel
  !   * Particle types
  !   * Xsecs
  !   * I/O Units
  !   * Physical constants
  !
  !=======================================================================================

  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  PUBLIC

  !--------------------------------------------------------- NUMERICAL CONSTANTS ----------
  ! Frequently used real numbers:
  real(real64), parameter ::  &
       &  zro              = 0.0_real64,                 &
       &  one              = 1.0_real64,                 &
       &  two              = 2.0_real64,                 &
       &  four             = 4.0_real64,                 &
       &  thousand         = 1000.0_real64
  !---------------------------------------------------------------------------------------


  !--------------------------------------------------------- PHYSICAL CONSTANTS ----------
  ! Particle masses.
  ! Ref. 1:  Current Coalescence value.
  ! Ref. 2:  Review of Particle Physics, Phys. Rev. D 66 (1 July 2002).
  ! Ref. 3:  CODATA Recommended Values of the Fundamental Physical Constants: 1998
  !          Peter J. Mohr and Barry N. Taylor
  !          National Institute of standards and Technology, Gaithersburg, MD 20899-8401
  !          http://physics.nist.gov/cuu/Constants/index.html
  !          For:  deuteron  triton  helion  alpha

  real(real64), parameter ::     &
       & neutron_mass  = 0.93958_real64,     &  != Ref. 1 ( 939.56533 MeV, Ref. 2)
       & proton_mass   = 0.938271998_real64, &  != Ref. 3 ( +/- 0.000038 MeV)
!               P.D.G. = 938.272013_real64
!                      ( 938.272 MeV, Ref. 2)
!                  2006: 938.27203 ( +/- 0.00008 )
       & avg_mass      = (neutron_mass + proton_mass)/two
  !---------------------------------------------------------------------------------------

end module coalescenceParams

