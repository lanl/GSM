
  function residualBaryons(residObj) result(numBaryons)

! ====================================================================
!
! Returns the number of baryons in a residObjual
!
! ====================================================================

    implicit none
    class(mDCMresidual), intent(in   ) :: residObj
    real(real64) :: numBaryons
    numBaryons = residObj%numBaryons
    return
! ====================================================================
  end function residualBaryons


  function residualProtons(residObj) result(numProtons)

! ====================================================================
!
! Returns the number of protons in a residual
!
! ====================================================================

    implicit none
    class(mDCMresidual), intent(in   ) :: residObj
    real(real64) :: numProtons
    numProtons = residObj%numProtons
    return
! ====================================================================
  end function residualProtons


  function residualKinEnergy(residObj) result(kinEnergy)

! ====================================================================
!
! Returns the kinetic energy (GeV) of a residual
!
! ====================================================================

    implicit none
    class(mDCMresidual), intent(in   ) :: residObj
    real(real64) :: kinEnergy
    kinEnergy = residObj%kinEnergy
    return
! ====================================================================
  end function residualKinEnergy


  function residualLinearMom(residObj) result(linearMom)

! ====================================================================
!
! Returns the linear moment of a residual
!
! ====================================================================

    implicit none
    class(mDCMresidual), intent(in   ) :: residObj
    real(real64), dimension(3) :: linearMom
    linearMom = residObj%linearMom
    return
! ====================================================================
  end function residualLinearMom


  function residualAngularMom(residObj) result(angularMom)

! ====================================================================
!
! Returns the angular momentum of a residual nucleus
!
! ====================================================================

    implicit none
    class(mDCMresidual), intent(in   ) :: residObj
    real(real64), dimension(3) :: angularMom
    angularMom = residObj%angularMom
    return
! ====================================================================
  end function residualAngularMom


  subroutine adjustresidual(residObj, numBaryons, numProtons, kinEnergy)

! ====================================================================
!
! This subroutine will adjust the internal characteristics of a
! residual nucleus object. This was created mostly to prevent
! accidental changes to it's values during development and for users.
!
! ====================================================================

    implicit none
    class(mDCMresidual), intent(inout) :: residObj      ! The residual object
    real(real64),        intent(in   ), optional :: numBaryons ! The new baryon number
    real(real64),        intent(in   ), optional :: numProtons ! The new proton number
    real(real64),        intent(in   ), optional :: kinEnergy  ! The new excitation energy [GeV]

    if (present(numBaryons))  residObj%numBaryons = numBaryons
    if (present(numProtons)) residObj%numProtons = numProtons
    if (present(kinEnergy)) residObj%kinEnergy = kinEnergy

    return
! ====================================================================
  end subroutine adjustresidual


  subroutine adjustMomentum(residObj, px, py, pz, x, y, z)

! ====================================================================
!
! This subroutine will adjust the internal momentum characteristics of a
! residual nucleus object. This was created mostly to prevent
! accidental changes to it's values during development and for users.
!
! ====================================================================

    implicit none
    class(mDCMresidual), intent(inout) :: residObj      ! The residObjual object
    real(real64),        intent(in   ) :: px         ! Linear momentum in X [units unknown]
    real(real64),        intent(in   ) :: py         ! Linear momentum in Y [units unknown]
    real(real64),        intent(in   ) :: pz         ! Linear momentum in Z [units unknown]
    real(real64),        intent(in   ) :: x          ! X position of the adjustment
    real(real64),        intent(in   ) :: y          ! Y position of the adjustment
    real(real64),        intent(in   ) :: z          ! Z position of the adjustment

    ! Adjust linear momentum
    residObj%linearMom(1) = residObj%linearMom(1) + px
    residObj%linearMom(2) = residObj%linearMom(2) + py
    residObj%linearMom(3) = residObj%linearMom(3) + pz

    ! Adjust angular momentum
    residObj%angularMom(1) = residObj%angularMom(1) + (z * py) - (y * pz)
    residObj%angularMom(2) = residObj%angularMom(2) + (x * pz) - (z * px)
    residObj%angularMom(3) = residObj%angularMom(3) + (y * px) - (x * py)

    return
! ====================================================================
  end subroutine adjustMomentum

