# ============================================================================ #
# GSM top-level CMakeLists.txt
#
# Options:
#   gsm.BUILD_DOCS    [ON]                Creates an API documentation build target
#   gsm.SIM_DIR       [ON]                Creates simulation directory build target
#   gsm.MPI           [OFF]               Builds GSM with MPI enabled
#   gsm.OpenMP        [OFF]               Builds GSM with OpenMP enabled
#   gsm.OpenACC       [OFF]               Builds GSM with OpenACC enabled
#   gsm.PACK          [on]                Packages GSM
#   gsm.SHARED        [BUILD_SHARED_LIBS] Builds GSM with shared libraries
#   gsm.TEST          [ON]                Enables testing of GSM via CTest
#
# Usage:
#    cmake [OPTIONS] ../
#
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)
project(gsm
  VERSION 1.1.0
  LANGUAGES Fortran
  )

# Set the CMake module directory
set(CMAKE_MODULE_PATH
  "${CMAKE_SOURCE_DIR}/.cmake"
  "${CMAKE_SOURCE_DIR}/.cmake/General"
  "${CMAKE_SOURCE_DIR}/.cmake/Utils"
  "${CMAKE_SOURCE_DIR}/.cmake/Compiler"
  "${CMAKE_SOURCE_DIR}/.cmake/CompilerCapabilities"
  "${CMAKE_SOURCE_DIR}/.cmake/FindLibrary"
  "${CMAKE_SOURCE_DIR}/.cmake/GSM"
  )

# Import modules for this CMakeLists
include (DefineProjectDefaults)
include (loadUtilities)
PreventInSourceBuild()

# Load project options
include (configureOptions)

# Obtain the compiler name and check its capabilities
include (configureCompiler)

# ============================================================================ #
#
# Configure source code
#
# ============================================================================ #
message(STATUS "${LightGreen}Configuring data and sub-models...${ColorReset}")
add_subdirectory(src)
add_subdirectory(Docs EXCLUDE_FROM_ALL)


# ============================================================================ #
#
# Configure executable
#
# ============================================================================ #
# NOTE: Simulation directory inclusion depends on this naming convention.
set (exec "x${PROJECT_NAME}${PROJECT_VERSION_MAJOR}")

# Set main executable's sources and dependencies
set(sources
  "${CMAKE_SOURCE_DIR}/src/main.f90"
  )
set(dependencies
  "GSM_Driver"
  )

# Add the exectable, link, and install it:
add_project_executable(${exec} "${sources}" "${dependencies}")

# ============================================================================ #

