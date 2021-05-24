! =============================================================================
!
!> \file
!> \brief  Contains the Loggers module
!> \author CMJ (XCP-3; LANL)
!
!> The Loggers module contains the Logger class used by GSM and its sub-models
!> and model-data. The base member variables, member-procedures, and constructors
!> are contained here.
!
! =============================================================================
module Loggers

    use, intrinsic:: iso_fortran_env, only: int32, output_unit
    use, intrinsic:: iso_C_binding,   only: c_int

    implicit none
    private

    public:: prepareMessage

    ! >>> PARAMETERIZED DEFAULTS
    include "parameters.f90"

    ! >>> OBJECT CONSTRUCTION
    include "constructorInterface.f90"

    ! >>> DATA OBJECTS
    include "logger.f90"

 contains

     ! >>> CONSTRUCTORS
     include "constructors.f90"

     ! >>> INTROSPECTION
     include "introspection.f90"

     ! >>> SETTERS
     include "setters.f90"

     ! >>> PROCEDURES
     include "procedures.f90"

end module Loggers

