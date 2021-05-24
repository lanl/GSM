# ============================================================================ #
#
# Sets the overarching compiler flags for the NAG compiler
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)

# Set preprocessor flag and flags common to the compiler
set (FPP_FLAG "-fpp")
set (common_flags "${FPP_FLAG}")

# Set the flags for various compilation types
set (CMAKE_Fortran_FLAGS_DEBUG   "\
${common_flags} -O0 \
-gline -C=all")
# workaround for nag 6.2
#set(CMAKE_Fortran_FLAGS_DEBUG "-C=array -C=alias -C=bits -C=calls -C=do -C=intovf -C=present -C=pointer -gline -O0")

set (CMAKE_Fortran_FLAGS_RELEASE "\
${common_flags} -O3")
