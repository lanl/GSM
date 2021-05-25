# ============================================================================ #
#
# Loads the configuration of the GSM packaging
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)

# Include package:
include(CPack)
message(${envDEBUG} "${Green}Configuring the ${PROJECT_NAME} package...${ColorReset}")


# Set basic package information
set(CPACK_PACKAGE_NAME "${PROJECT_NAME}${VERSION_MAJOR}-${VERSION_MINOR}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
  "The Generalized Spallation Model Event Generator")
set(CPACK_PACKAGE_VENDOR
  "Los Alamos National Laboratory and Idaho State University")
set(CPACK_RESOURCE_FILE "${CMAKE_SOURCE_DIR}/ReadMe.md")


# Set other package information:
set(CPACK_PACKAGE_DIRECTORY "${CMAKE_BINARY_DIR}/package")
set(CPACK_STRIP_FILES TRUE)

# ============================================================================ #

