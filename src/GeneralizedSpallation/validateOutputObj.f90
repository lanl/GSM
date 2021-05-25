
! ======================================================================
!
! Validates the GSMOutput object
!
!
! Writteny by CMJ, XCP-3 (07/2019)
!
! ======================================================================

  subroutine validateOutputObj(output)

    use, intrinsic:: iso_fortran_env, only: int32

    implicit none
    class(GSMOutput), intent(inout) :: output

    integer(int32) :: iter = 0_int32

! ======================================================================

    ! Validate num. inelastic events
    if(output%numInelasticEvents < 0) output%numInelasticEvents = &
         & defaultInelasticEvents

    ! Validate what to calculate
    if(output%angularSpectra < 0 .or. output%angularSpectra > 2) &
         & output%angularSpectra = defaultAngSpec
    if(output%energySpectra < 0 .or. output%energySpectra > 2) &
         & output%energySpectra = defaultEnergySpec
    if(output%doubleDiffSpectra < 0 .or. output%doubleDiffSpectra > 2) &
         & output%doubleDiffSpectra = defaultDoubDiffSpec
    if(output%multiplicities < 0 .or. output%multiplicities > 2) &
         output%multiplicities = defaultMultiplicity
    if(output%channelCrossSection < 0 .or. output%channelCrossSection > 1) &
         & output%channelCrossSection = defaultChannelXSCalc
    if(output%nuclideCrossSection < 0 .or. output%nuclideCrossSection > 3) &
         & output%nuclideCrossSection = defaultNuclideXSCalc

    ! Validate angle bin widths
    if(output%deltaTheta < 0) output%deltaTheta = defaultDTheta
    if(output%pisaDTheta < 0) output%pisaDTheta = defaultDTheta

    ! Validate angle bins
    angles: do iter = 1_int32, numAnglesConsidered
       ! PISA:
       if(output%pisaAngles(iter) < 0 .or. output%pisaAngles(iter) >= 360) &
            & output%pisaAngles(iter) = defaultPisaAngles(iter)

       ! Normal bins:
       if(output%angleBins(iter)%lowerBound < 0 .or. &
            & output%angleBins(iter)%lowerBound >= 360) &
            & output%angleBins(iter)%lowerBound = defaultAngleBins(iter)%lowerBound
       if(output%angleBins(iter)%upperBound < 0 .or. &
            & output%angleBins(iter)%upperBound >= 360) &
            & output%angleBins(iter)%upperBound = defaultAngleBins(iter)%upperBound
       if(output%angleBins(iter)%upperBound <= output%angleBins(iter)%lowerBound) then
          output%angleBins(iter)%lowerBound = defaultAngleBins(iter)%lowerBound
          output%angleBins(iter)%upperBound = defaultAngleBins(iter)%upperBound
       end if
    end do angles

    energy: do iter = 1_int32, numEnergyBins
       ! Sub-step
       if(output%energyBinSubStep(iter) < 0.0_real64) &
            & output%energyBinSubStep(iter) = defaultEnergyBinSubStep(iter)

       ! Bins
       if(output%energyBins(iter)%lowerBound < 0.0_real64 ) &
            & output%energyBins(iter)%lowerBound = defaultEnergyBins(iter)%lowerBound
       if(output%energyBins(iter)%upperBound < 0.0_real64 ) &
            & output%energyBins(iter)%upperBound = defaultEnergyBins(iter)%upperBound
       if(output%energyBins(iter)%upperBound <= output%energyBins(iter)%lowerBound) then
          output%energyBins(iter)%lowerBound = defaultEnergyBins(iter)%lowerBound
          output%energyBins(iter)%upperBound = defaultEnergyBins(iter)%upperBound
       end if
    end do energy

    return
! ======================================================================
  end subroutine validateOutputObj
