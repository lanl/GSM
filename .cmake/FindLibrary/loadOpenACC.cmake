# ============================================================================ #
#
# Finds and sets options for OpenACC usage
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.10.0)

if (${PROJECT_NAME}.OpenACC)
  find_package(OpenACC)

  if (OPENACC_FOUND OR
      (${CMAKE_VERSION} VERSION_LESS "3.25" AND
      NOT "${OpenACC_Fortran_VERSION}" STREQUAL ""))
    # Print that OpenACC was found (similar to how we would for higher CMake
    # versions)
    if (${CMAKE_VERSION} VERSION_LESS "3.25")
      message("-- Found OpenACC: TRUE (found version \"${OpenACC_Fortran_VERSION}\")")
    endif()

    # CMake 3.16+ defines the compiler flags; they may not exist for
    # lesser versions. We'll assume them since it's pretty consistent
    if (${CMAKE_VERSION} VERSION_LESS "3.16")
      set (OpenACC_Fortran_FLAGS "-fopenacc")
    endif()

    # Add OpenACC flags to debug/release flags now
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
