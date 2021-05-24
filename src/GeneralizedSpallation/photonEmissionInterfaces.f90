
  abstract interface 
     ! PREEQUILIBRIUM PHOTON EMISSION INTERFACE:
     subroutine preequilibriumPHOTONEMISSION(preeqFrag, residual)
       use preequilibriumClass, only: preequilibriumFragment, residualNucleus
       implicit none
       type(preequilibriumFragment), intent(in) :: preeqFrag
       type(residualNucleus),        intent(in) :: residual
     end subroutine preequilibriumPHOTONEMISSION
     ! EVAPORATION PHOTON EMISSION INTERFACE:
     subroutine evaporationPHOTONEMISSION(evapFrag, residA, residZ, residKE)
       use, intrinsic:: iso_fortran_env, only: int32, real64
       use evaporationClass, only: evaporationFragment
       implicit none
       type(evaporationFragment), intent(in) :: evapFrag
       real(real64),              intent(in) :: residA
       real(real64),              intent(in) :: residZ
       real(real64),              intent(in) :: residKE
     end subroutine evaporationPHOTONEMISSION
     ! GSM PHOTON EMISSION INTERFACE:
     subroutine PHOTOEMISSION(fragment, residual)
       import GSMProgeny, GSMResidual
       implicit none
       type(GSMProgeny), intent(in) :: fragment
       type(GSMResidual), intent(in) :: residual
     end subroutine PHOTOEMISSION
  end interface
