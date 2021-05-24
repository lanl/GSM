# ============================================================================ #
#
# Sets the overarching compiler flags for the Intel compiler
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)

# Set preprocessor flag and flags common to the compiler
set (FPP_FLAG "-fpp")
set (common_flags "${FPP_FLAG} \
-free -stand f08")

# Set compiler flags for each compilation type
set (CMAKE_Fortran_FLAGS_RELEASE "\
${common_flags} -O3")

set (CMAKE_Fortran_FLAGS_DEBUG   "\
${common_flags} -O0 \
-traceback \
-check uninit -save-temps")

# ============================================================================ #

