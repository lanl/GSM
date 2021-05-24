# ============================================================================ #
#
# Configures the various build variables as appropriate
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)

# Set locations for all built files:
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Set messages for installation
set(CMAKE_INSTALL_MESSAGE ALWAYS)

