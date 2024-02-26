# ============================================================================ #
#
# Configures the compilation of the project by setting flags and testing the
# compiler's features.
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)
PreventInSourceBuild()

# Obtain the compiler name and check capabilities
include(${CMAKE_Fortran_COMPILER_ID} RESULT_VARIABLE found)
if(NOT found)
  message(WARNING
    "An unrecognized Fortran compiler is being used.\n"
    "   ${PROJECT_NAME} may not compile successfully.")
  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -cpp")
  set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -cpp")
endif()

# Print the flags being used to the user

message(${envDEBUG}
  "${PROJECT_NAME} will be compiled with the following flag sets:\n"
  "     (Release) ${Cyan}${CMAKE_Fortran_FLAGS_RELEASE}${ColorReset}\n"
  "     (Debug)   ${Cyan}${CMAKE_Fortran_FLAGS_DEBUG}${ColorReset}")

# Check if the compiler has the capabilities required
include (CheckCompilerCapabilities)

# ============================================================================ #
