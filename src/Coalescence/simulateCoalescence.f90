
  subroutine simulateCoalescence( coalObj, results )

! ==============================================================================
!
! This subroutine is called by client programs to simulate coalescence physics.
! This subroutine acts as a "layer" to the primary routine, 'coalesl'.
!
! The primary purpose of this layer is to create several integers that can be
! used to create better memory management and ensure that arrays internal to
! 'coalesl' do not get exceeded.
!
!
! Written by CMJ, XCP-3, 01/2019
!
! ==============================================================================

    use, intrinsic:: iso_fortran_env, only: int32

    implicit none
    class(Coalescence),       intent(inout) :: coalObj
    type(coalescenceResults), intent(inout) :: results

    ! Integers for the number of fragments that can form w/ A= [2, 3, 4, 6, 7]
    integer(int32) :: numA2, numA3, numA4, numA6, numA7

! ==============================================================================

    ! Ensure class was constructed:
    if ( .not. coalObj%constructed ) then
       ! Class hasn't been constructed; print error and stop                                   
       write(coalObj%io%message, 3000)
       call coalObj%io%print(1, 2, coalObj%io%message)
       write(coalObj%io%message, 3999)
       call coalObj%io%print(1, 2, coalObj%io%message)
       results%simState = 10   ! Flag that simulation did not complete due to error
       return
    end if


    ! Check that value of numFragments is valid
    if ( results%numParticles > results%partBnkSize ) then
       ! Warning - there cannot be more particles than allowed in the arrays. Warn client
       write(coalobj%io%message, 3100) results%numParticles, results%partBnkSize
       call coalObj%io%print(2, 3, coalObj%io%message)
       results%numParticles = results%partBnkSize
       write(coalobj%io%message, 3110) results%partBnkSize
       call coalObj%io%print(2, 3, coalObj%io%message)
       results%simState = 1
    else if ( results%numParticles < 0 ) then
       ! Warning - an unphysical number of particles was requested; correct.
       write(coalobj%io%message, 3150) results%numParticles, results%partBnkSize
       call coalObj%io%print(2, 3, coalObj%io%message)
       results%numParticles = results%partBnkSize
       results%simState = 2
    end if


    ! Ensure results object was constructed
    if ( .not. results%constructed ) then
       ! Results object wasn't constructed; print error and stop                              
       write(coalObj%io%message, 3200)
       call coalObj%io%print(1, 2, coalObj%io%message)
       write(coalObj%io%message, 3999)
       call coalObj%io%print(1, 2, coalObj%io%message)
       results%simState = 11   ! Flag that simulation did not complete due to error
       return
    end if


    ! Check if anything can coalesce:
    if ( results%numParticles <= 1 ) return


    ! Determine maximum sizes for various compound fragment arrays
    numA2 = results%numParticles / 2
    numA3 = results%numParticles / 3
    numA4 = results%numParticles / 4
    numA6 = results%numParticles / 6
    numA7 = results%numParticles / 7


    ! Perform coalescence simulation
    call coalObj%coalesl( results, numA2, numA3, numA4, numA6, numA7 )


    return
! ==============================================================================
3000 format("The coalescence object has not been ", &
          & "constructed.")
3100 format("Cannot coalesce ", i6, " particle(s). No more than ", i6, &
          & " particle(s) are permitted in memory.")
3110 format("   Limiting to ", i6, " particles.")
3150 format("Cannot coalesce ", i6, " particle(s). Assuming ", i6, &
          & " particle(s) exist instead.")
3200 format("The results object passed in to the coalescence object ", &
          & "was not constructed.")
3999 format("   Unable to simulate coalescence physics.")
! ==============================================================================
  end subroutine simulateCoalescence
