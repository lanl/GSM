
  function coalesResultsConstructorMain( clientParticleBank, numParticles ) &
       & result(results)

! ==============================================================================
!
! This function is used to construct a results object for clients to use.
! Clients simply pass in the particle bank (using the coalescence particle type)
! and the number of fragments the client wishes to consider for coalescence
!
!
! Written by CMJ, XCP-3, July 2018
! Edited  by CMJ, XCP-3, Jan. 2019
!
! ==============================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    implicit none

    type(coalescenceParticle), intent(in), dimension(:), target   :: clientParticleBank
    integer(int32),            intent(in),               optional :: numParticles
    type(coalescenceResults)                                      :: results

! ==============================================================================

    ! Set address of secondary arrays
    results%partBnk => clientParticleBank


    ! Set amount of particles
    results%partBnkSize  = size( results%partBnk )
    if ( present(numParticles) ) then
       ! NOTE: This is checked for errors prior to simulation
       results%numParticles = numParticles
    else
       results%numParticles = results%partBnkSize
    end if


    ! The results object was constructed successfully
    results%constructed = .TRUE.

    return
! ==============================================================================
  end function coalesResultsConstructorMain
