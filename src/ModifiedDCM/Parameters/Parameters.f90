! copyright LANS/LANL/DOE - see file COPYRIGHT_INFO

module modifiedDCMParams
  !=======================================================================================
  ! Description:
  ! Various parameters used by many routines in modifiedDCMParams.
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

  use, intrinsic :: iso_fortran_env, only: int32, real64
  use numbers
  implicit none
  PUBLIC

  !--------------------------------------------------------- PHYSICAL CONSTANTS ----------
  ! Values cited as PDG 2011 are from the following reference:
  !   K. Nakamura et al. (Particle Data Group), Journal of Physics G37, 075021 (2010)
  !      and 2011 partial update for the 2012 edition.
  ! Fundamental physical constants:
  real(real64), parameter ::  &
       & fscon         = 137.0393_real64,             & != Inverse fine-structure constant (dimensionless).
       & fscon_inv     = 137.035999679_real64,        & != Inverse fine-structure constant (dimensionless). PDG 2011.
       & hbarc         = 197.327053_real64              != hbar * c or hbar/c?


  ! Particle masses.
  ! Ref. 1:  Current modifiedDCMParams value.
  ! Ref. 2:  Review of Particle Physics, Phys. Rev. D 66 (1 July 2002).
  ! Ref. 3:  CODATA Recommended Values of the Fundamental Physical Constants: 1998
  !          Peter J. Mohr and Barry N. Taylor
  !          National Institute of standards and Technology, Gaithersburg, MD 20899-8401
  !          http://physics.nist.gov/cuu/Constants/index.html
  !          For:  deuteron  triton  helion  alpha

  real(real64), parameter ::     &
       & electron_mass = 0.511008_real64,   &  != Current modifiedDCMParams value. However, Review of Particle Physics,
!               P.D.G.OB = 0.51099891_real64
                                          != Phys. Rev. D 66 (1 July 2002) gives 0.5109989 MeV.
       & neutron_mass  = 939.58_real64,     &  != Ref. 1 ( 939.56533 MeV, Ref. 2)
       & proton_mass   = 938.271998_real64, &  != Ref. 3 ( +/- 0.000038 MeV)
!               P.D.G. = 938.272013_real64
!                      ( 938.272 MeV, Ref. 2)
!                  2006: 938.27203 ( +/- 0.00008 )
       & amu           = 931.4943d0,      & ! Average amu mass (value in CEM)
       & nucleon_mass  = 931.4943_real64,   & ! Average amu mass
       & nucleon_mass_estimate = 940.0_real64, & ! poor estimation that is used in calculation
       & aneut         = 1.008664967_real64,   & != Neutron mass in a.m.u.
       & emnucm        = 938.919d0,          & != Nucleon mass [MeV/c^2]
       & emnucg        = 0.938919d0,         & != Nucleons mass [GeV/c^2]
       & emnucb        = 0.9315014d0,        & != Average nucleon mass [GeV/c^2]
       & emnuct        = 0.9314943d0,        & != Nucleon mass
       & emneut        = 0.9395656d0,        & != Neutron mass
       & emprot        = 0.9382723d0,        & != Proton mass
       & empich        = 0.139568d0,         & != Pion mass
       & empi0         = 0.134973d0            != Pi0 mass


  ! Derived physical constants:
  real(real64), parameter ::  &
       & alpha_fsc       = one / fscon       != Fine-structure constant (dimensionless).
  !---------------------------------------------------------------------------------------

end module modifiedDCMParams

