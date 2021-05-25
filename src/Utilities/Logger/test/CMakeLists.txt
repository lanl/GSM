# ============================================================================ #
#
# CMakeLists.txt for testing the Logger object
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)
PreventInSourceBuild()

# Set sub-model's name
set(test_name "Logger")

# ============================================================================ #

# Load the sources used for the test
set(sources
  ${CMAKE_CURRENT_LIST_DIR}/main.f90
  ${CMAKE_CURRENT_LIST_DIR}/tstLogger.f90
  )

# List the test's dependencies
set(dependencies
  Logger
  )

# ============================================================================ #

# Create the test's executable
add_test_executable(${test_name} "${sources}" "${dependencies}")

# ============================================================================ #
