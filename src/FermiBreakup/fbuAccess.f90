
! ====================================================================
!
! This file contains the various access methods by which clients may
! query the internal state of the FBU instance.
!
!
! Written by CMJ, XCP-3 (03/2019)
!
! ====================================================================

  function recommendFBUReal(fbuObj, numBaryonsReal) result(simFBU)

! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    class(FermiBreakUp), intent(in   ) :: fbuObj
    real(real64),        intent(in   ) :: numBaryonsReal
    logical                            :: simFBU

    integer(int32) :: numBaryonsInt

! ====================================================================

    numBaryonsInt = nint(numBaryonsReal)
    simFBU = fbuObj%recommendFBUInt(numBaryonsInt)

    return
! ====================================================================
  end function recommendFBUReal


  function recommendFBUInt(fbuObj, numBaryons) result(simFBU)

! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32

    implicit none
    class(FermiBreakUp), intent(in   ) :: fbuObj
    integer(int32),      intent(in   ) :: numBaryons
    logical                            :: simFBU

! ====================================================================

    simFBU = (numBaryons <= fbuObj%options%recNumNucleons ) .and. &
         &   (numBaryons >= minAllowedNucleons)

    return
! ====================================================================
  end function recommendFBUInt


  function properlyConstructed(fbuObj) result(constructed)

! ====================================================================
!
! Returns whether or not the object is constructed
!
! ====================================================================

    implicit none
    class(FermiBreakUp), intent(in   ) :: fbuObj
    logical :: constructed

! ====================================================================

    constructed = fbuObj%constructed

    return
! ====================================================================
  end function properlyConstructed


  function queryOptions(fbuObj) result(options)

! ====================================================================
!
! Procedure returns the options object utilized by the FBU object
!
! ====================================================================

    implicit none
    class(FermiBreakUp), intent(in   ) :: fbuObj
    type(fermiBreakUpOptions) :: options

! ====================================================================

    options = fbuObj%options

    return
! ====================================================================
  end function queryOptions
