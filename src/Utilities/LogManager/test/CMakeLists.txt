# ============================================================================ #
#
# CMakeLists.txt for testing the LogManager object
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)
PreventInSourceBuild()

# Set sub-model's name
set(test_name "LogManager")

# ============================================================================ #

# Load the sources used for the sub-model
set(sources
  ${CMAKE_CURRENT_LIST_DIR}/main.f90
  ${CMAKE_CURRENT_LIST_DIR}/tstLogManager.f90
  )

# List the sub-model's dependencies
set(dependencies
  Logger
  LogManager
  )

# ============================================================================ #

# Create the sub-model's library
add_test_executable(${test_name} "${sources}" "${dependencies}")

# ============================================================================ #
