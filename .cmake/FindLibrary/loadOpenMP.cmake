# ============================================================================ #
#
# Finds and sets options for OpenMP usage
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)

if (${PROJECT_NAME}.OpenMP)
  find_package(OpenMP)
  if (OPENMP_FOUND)
    # Append OpenMP compiler flags:
    set(CMAKE_Fortran_FLAGS_DEBUG
      "${CMAKE_Fortran_FLAGS_DEBUG} ${OpenMP_Fortran_FLAGS}")
    set(CMAKE_Fortran_FLAGS_RELEASE
      "${CMAKE_Fortran_FLAGS_RELEASE} ${OpenMP_Fortran_FLAGS}")
  else ()
    set(${PROJECT_NAME}.OpenMP OFF)
    message(${envNOTICE} "OpenMP suport not available")
  endif()
endif ()


# ============================================================================ #
