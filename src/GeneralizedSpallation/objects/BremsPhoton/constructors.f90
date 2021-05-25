! =============================================================================
!
!> \file
!> \brief  Contains the BremsPhoton object constructors
!> \author CMJ (XCP-3; LANL)
!
! =============================================================================

! =============================================================================
!
!> \fn    constructorMain
!> \brief Main Logger constructor
!
!> Main constructor for the BremsPhoton object. Clients may use this to
!> contain characteristics of a brems. photon, specifying an energy range to
!> sample from
!
! ARGUMENTS:
!> \param[in   ] minE        Minimum energy of the brems. photon
!> \param[in   ] maxE        Maximum energy of the brems. photon
!
! =============================================================================
 function constructorMain( &
         &  minE &
         & ,maxE &
         & ) &
         & result(this)

     use, intrinsic:: iso_fortran_env, only: real64

     implicit none
     real(real64), intent(in   ):: minE
     real(real64), intent(in   ):: maxE
     type(BremsPhoton):: this

! =============================================================================

     ! Interface to establish internal values of the object
     call this%setTMin(minE)
     call this%setTMax(maxE)
     call this%resetTEqv()
     call this%resetSXAbs()

     return
! =============================================================================
end function constructorMain


