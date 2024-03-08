# ============================================================================ #
#
# Finds and sets options for OpenACC usage
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.10.0)

if (${PROJECT_NAME}.OpenACC)
  find_package(OpenACC)
  if (OPENACC_FOUND)
    # Append OpenACC compiler flags:
    set(CMAKE_Fortran_FLAGS_DEBUG
      "${CMAKE_Fortran_FLAGS_DEBUG} ${OpenACC_Fortran_FLAGS}")
    set(CMAKE_Fortran_FLAGS_RELEASE
      "${CMAKE_Fortran_FLAGS_RELEASE} ${OpenACC_Fortran_FLAGS}")
  else ()
    set(${PROJECT_NAME}.OpenACC OFF)
    message(${envNOTICE} "OpenACC suport not available")
  endif()
endif ()


# ============================================================================ #
