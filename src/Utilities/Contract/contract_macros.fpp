! =============================================================================
!
!> \file
!> \brief  Contains the macros used to simplify calling various Contracts
!> \author CMJ (XCP-3; LANL)
!
! Defines various macros to simplify calling and inlining for performance
! Simply "#include" this in the consuming module, and any sub-files using these
! will also need "#include"d within the module.
!
! =============================================================================

#ifndef contract_macros_f90
#define contract_macros_f90

! Define various contract levels
#define PRODUCTION_LEVEL    0
#define PRECONDITION_LEVEL  1
#define INTRASCOPE_LEVEL    2
#define POSTCONDITION_LEVEL 3
#define DEBUG_LEVEL         4

! Define the Contracts level for the build
#ifndef CONTRACTS_LEVEL
#define CONTRACTS_LEVEL DEBUG_LEVEL
#endif


! Define "Require" - used for parameter validation (pre-condition checks)
#if CONTRACTS_LEVEL >= 1
#define Require(COND) call require(COND, __FILE__, __LINE__)
#else
#define Require(COND)
#endif


! Define "Check" - used for intra-scope validation (mid-calculation checks)
#if CONTRACTS_LEVEL >= 2
#define Check(COND) call check(COND, __FILE__, __LINE__)
#else
#define Check(COND)
#endif


! Define "Ensure" - used for post-condition validation (after-calculation checks)
#if CONTRACTS_LEVEL >= 3
#define Ensure(COND) call ensure(COND, __FILE__, __LINE__)
#else
#define Ensure(COND)
#endif


! Define "Insist" - always on for conditions that must always be met.
#define Insist(COND, MSG) \
   call insist(COND, MSG, __FILE__, __LINE__) \


#define NotImplemented(FEATURE) \
   call not_implemented(FEATURE, __FILE__, __LINE__)
 

#endif
