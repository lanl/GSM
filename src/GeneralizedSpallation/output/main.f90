
! =============================================================================
!
!> \file
!> \brief   Contains the main include files for all output-file related
!>          procedures
!> \author  CMJ, XCP-3 (LANL)
!
! ==============================================================================

   !> Writes the license and reaction information to the output file
   include "output/prinp.f90"

   !> Main output generating procedure
   include "output/typeout.f90"

   !> Estimates inverse cross section for simulation
   include "output/siginv.f90"

   !> Other procedures that write distribution outputs
   include "output/pdisnm.f90"
   include "output/propan.f90"
   include "output/prrdis.f90"
   include "output/prtdadz.f90"
   include "output/prtdist.f90"
   include "output/prtmult.f90"

