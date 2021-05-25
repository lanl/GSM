# ============================================================================ #
#
# Checks if the Fortran compiler supports features utilized by the project
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)

# Include the function to test features
include (CheckFortranSource)


# ============================================================================ #
# Test various compiler features now:
# ============================================================================ #
Check_Fortran_Source_Compile (
  ${CMAKE_CURRENT_LIST_DIR}/supportForISOFortranEnv.F90
  SUPPORT_FOR_ISO_FORTRAN_ENV
)


Check_Fortran_Source_Compile (
  ${CMAKE_CURRENT_LIST_DIR}/supportForISOCBinding.F90
  SUPPORT_FOR_ISO_C_BINDING
)


Check_Fortran_Source_Compile (
  ${CMAKE_CURRENT_LIST_DIR}/supportForAssumedType.F90
  SUPPORT_FOR_ASSUMED_TYPE
)


Check_Fortran_Source_Compile (
  ${CMAKE_CURRENT_LIST_DIR}/supportForDynamicStrings.F90
  SUPPORT_FOR_DYNAMIC_STRINGS
)


Check_Fortran_Source_Compile (
  ${CMAKE_CURRENT_LIST_DIR}/supportForPreprocessor.F90
  SUPPORT_FOR_FORTRAN_PREPROCESSOR
)


# Check_Fortran_Source_Compile (
#   ${CMAKE_CURRENT_LIST_DIR}/supportForSelectType.F90
#   SUPPORT_FOR_FORTRAN_SELECT_TYPE
# )



