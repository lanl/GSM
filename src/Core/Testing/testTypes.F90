! =============================================================================
!>
!> \file
!> \brief  Contains the Testing types
!> \author CMJ
!>
! =============================================================================

  !> \brief A basic test object
  !>
  !> A basic test object. It provides a subroutine for the test and indicates
  !> a pass or fail state of the test.
  type, public :: Test
     private

     !> \brief Description of the test
     character(:), private, allocatable :: description

     !> \brief Status of the test
     type(TestStatus), private :: status = UNTESTED

     !> \brief Procedure to execute for the test


   contains
     private

     !> \brief Executes the test
     procedure, public :: execute
     generic,   public :: run => execute

     !> \brief Returns if the test passed or not
     procedure, public :: passed
     generic,   public :: succeeded => passed
     procedure, public :: failed
     generic,   public :: failure => failed

     !> \brief Returns the description of the test
     procedure, public :: description
     generic,   public :: describe => description

     !> \brief Sets of the description of the test
     procedure, public :: setDescription

     !> \brief Sets the procedure to execute for testing
     procedure, public :: setTestProcedure
     generic,   public :: setProcedure

  end type Test


  !> \brief A set of basic test objects.
  !>
  !> The most basic testing object. Real tests should be derived from this
  !> object as it provides the basic framework for testing.
  type, public, abstract :: TestGroup
     private

     !> \brief Brief description of the test set
     character(:), private, allocatable :: description
     
     !> \brief Set of sub-tests contained by the parent test
     class(Test), private, dimension(:), allocatable :: tests

   contains
     private
     
  end type TestGroup
       
