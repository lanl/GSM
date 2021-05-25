! =============================================================================
!
!> \file
!> \brief  Contains the implementation of the Attribute object
!> \author CMJ (XCP-3; LANL)
!
! =============================================================================
!
! DATA TYPE DESCRIPTION
!
!> \class Attribute
!> \brief Defines the base Attribute object
!
!> Defines the base Attribute object used by GSM and its sub-models and
!> model-data.
!
! =============================================================================
type, public:: Attribute
    private
    
    !> \brief Flags the construction state of the object
    logical, private:: b_constructed = defaultConstructed

    !> \brief Name of the attribute object
    character, private, allocatable:: b_name

    !> \brief Description of the attribute object
    character, private, allocatable:: b_description

    !> \brief Message logging object
    type(Logger), private, pointer:: b_log => NULL()

 contains
    private

    ! >>> INTROSPECTION

    !> \brief Returns if the object was constructed or not
    procedure, public:: getConstructed

    !> \brief Returns the name of the object
    procedure, public:: getName

    !> \brief Returns the description of the object
    procedure, public:: getDescription

    !> \brief Returns the object's Logger
    procedure, public:: getLogger


    ! >>> SETTERS

    !> \brief Sets the construction state of the object
    procedure, private:: setConstructed

    !> \brief Sets the name of the object
    procedure, public:: setName

    !> \brief Sets the object's description
    procedure, public:: setDescription

    !> \brief Sets the Logger contained by the object
    procedure, public:: setLogger


    ! >>> GENERIC PROCEDURES

    !> \brief Returns or sets the name of the object
    generic, public:: name => getName

    !> \brief Returns or sets the description of the object
    generic, public:: description => getDescription

    !> \brief Returns or sets the Logger object
    generic, public:: log => getLogger

    !> \brief Returns the construction flag of the object
    generic, public:: constructed => getConstructed

end type Attribute

