# ============================================================================ #
#
# Finds and sets options for MPI usage
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.10.0)

if (${PROJECT_NAME}.MPI)
  find_package(MPI)
  if (MPI_FOUND)
    include_directories (${MPI_Fortran_INCLUDE_DIRS})
  else ()
    set(${PROJECT_NAME}.MPI OFF)
    message(${envNOTICE} "MPI suport not available")
  endif()
endif ()


# ============================================================================ #
