# ============================================================================ #
#
# Loads the configuration of the GSM code coverage
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)

# Load code coverage module IF a debug build
if (CMAKE_BUILD_TYPE AND
    CMAKE_BUILD_TYPE MATCHES "^(Debug|DEBUG)$")
  message(${envDEBUG} "${PROJECT_NAME} code coverage analysis is enabled.")
  include (CodeCoverage)
  append_coverage_compiler_flags()
else()
  # For non-debug builds, simply don't allow coverage as an option
  message(WARNING
    "Code Coverage is not allowed for non-DEBUG builds.")
endif()

# ============================================================================ #
#
