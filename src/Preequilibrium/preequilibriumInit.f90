
  subroutine preequilibriumInit()

! ======================================================================
!
! Sets up preequilibrium data
!
! ======================================================================

    use preequilibriumData, only: initPreequilibriumData, preeqDataDeclared
    implicit none

! ======================================================================

    ! Set data used by the preequilibrium model if not yet done
    if ( .not. preeqDataDeclared ) call initPreequilibriumData()
    return

! ======================================================================
  end subroutine preequilibriumInit
