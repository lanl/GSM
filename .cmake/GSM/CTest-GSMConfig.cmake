# ============================================================================ #
#
# Loads the configuration of the GSM testing
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)

# Load testing
include (CTest)
# enable_testing()


add_custom_target(build_and_test
  ${CMAKE_CTEST_COMMAND}
  COMMENT "${Cyan}Building and running tests...${ColorReset}"
  )
message(${envDEBUG} "${PROJECT_NAME} testing is enabled.")
if (NOT "${CMAKE_BUILD_TYPE}" MATCHES "^(Debug|DEBUG)$")
  message(${envDEBUG}
    "   Set CMAKE_BUILD_TYPE to DEBUG when running tests for a robust usage.")
endif ()


# Set basic testing information
set(CTEST_PROJECT_NAME "${PROJECT_NAME}}")
set(CTEST_NIGHTLY_START_TIME "00:00:00 MDT")
set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "UNSET")
set(CTEST_DROP_LOCATION "UNSET")
set(CTEST_DROP_SITE_CDASH TRUE)

# ============================================================================ #
#
