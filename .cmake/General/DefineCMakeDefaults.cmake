# ============================================================================ #
#
# Defines some default CMake values
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)


# Always include srcdir and builddir in include path
# This saves typing ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY} in
# about every subdir
# since cmake 2.4.0
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Put the include dirs which are in the source or build tree
# before all other include dirs, so the headers in the sources
# are prefered over the already installed ones
# since cmake 2.4.1
set(CMAKE_INCLUDE_DIRECTORIES_PROJECT_BEFORE ON)

# Use colored output
# since cmake 2.4.0
set(CMAKE_COLOR_MAKEFILE ON)

# Set the default build type to release with release info
set(PROJECT_CONFIGURATIONS
  DEBUG
  RELEASE
  )
set(SUPPORTED_CONFIGURATIONS "${PROJECT_CONFIGURATIONS}" CACHE STRING
  "Supported ${PROJECT_NAME} build type configurations.")

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RELEASE" CACHE STRING "Choose the type of build." FORCE)
  message(STATUS "Setting build type to '${CMAKE_BUILD_TYPE}' as none was specified.")
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "DEBUG" "RELEASE")
endif()
