# ============================================================================ #
#
# Sets the overarching compiler flags for the GNU compiler
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)

# Set preprocessor flag and flags common to the compiler
set (FPP_FLAG "-cpp")
set (common_flags "\
${FPP_FLAG} \
-fmodule-private \
-fimplicit-none \
-ffree-line-length-128")

# Set the flags for various compilation types
set (CMAKE_Fortran_FLAGS_RELEASE "\
${CMAKE_Fortran_FLAGS_RELEASE} ${common_flags} -O3")

set (CMAKE_Fortran_FLAGS_DEBUG "\
${CMAKE_Fortran_FLAGS_DEBUG} \
${common_flags} -O0 \
-Wall -Wno-conversion \
-Wno-line-truncation -Wno-error=line-truncation \
-pedantic -fcheck=all -fbacktrace \
-finit-real=snan")
if (NOT "${CMAKE_Fortran_COMPILER_VERSION}" VERSION_LESS "5")
  set (CMAKE_Fortran_FLAGS_DEBUG "\
${CMAKE_Fortran_FLAGS_DEBUG} \
-ffpe-summary=all")
endif()

# set (CMAKE_Fortran_FLAGS_DEBUG "\
# ${common_flags} -O0 \
# -Wall -Wno-conversion \
# -pedantic -fcheck=all -fbacktrace \
# -ffpe-trap=invalid,zero,overflow \
# -finit-real=snan")

# ============================================================================ #

