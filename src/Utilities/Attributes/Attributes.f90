! =============================================================================
!
!> \file
!> \brief  Contains the Attributes module
!> \author CMJ (XCP-3; LANL)
!
!> The Attributes module contains the Attribute class used by GSM and its
!> sub-models and model-data. The base member-variables, member-procedures, and
!> constructors are contained here.
!
! =============================================================================
module Attributes

    use Loggers, only: Logger

    implicit none
    private

    ! >>> PARAMETERIZED DEFAULTS
    include "parameters.f90"


    ! >>> OBJECT CONSTRUCTION
    include "constructorInterface.f90"


    ! >>> DATA OBJECTS
    include "attribute.f90"

 contains

     ! >>> CONSTRUCTORS
     include "constructors.f90"

     ! >>> INTROSPECTION
     include "introspection.f90"

     ! >>> SETTERS
     include "setters.f90"

end module Attributes

