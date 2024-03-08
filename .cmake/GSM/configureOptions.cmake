# ============================================================================ #
#
# Loads the options available to GSM
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)

# Define the function to add an option to the GSM project
function(add_option
    option_name description default)
  set(opt_name "${PROJECT_NAME}.${option_name}")
  option(${opt_name} "${description}" ${default})
endfunction()

# ============================================================================ #

# Now specify the available options:
add_option("BUILD_DOCS"
  "Create and install the HTML based API documentation for ${PROJECT_NAME} \
  (requires Doxygen)"
  ON)
add_option("MPI"
  "Compile ${PROJECT_NAME} with MPI"
  OFF)
add_option("OpenMP"
  "Compile ${PROJECT_NAME} with OpenMP"
  OFF)
add_option("OpenACC"
  "Compile ${PROJECT_NAME} with OpenACC"
  OFF)
add_option("PACKAGE"
  "Package ${PROJECT_NAME}"
  ON)
add_option("SHARED"
  "Build shared libraries"
  "${BUILD_SHARED_LIBS}")
add_option("TEST"
  "Enable the unit testing in GSM"
  ON)


# \todo Deprecate this option
add_option("SIM_DIR"
  "Make a simulation directory"
  ON)

# ============================================================================ #
# Perform any special task for the options
# ============================================================================ #

if(${PROJECT_NAME}.BUILD_DOCS)
  find_package(Doxygen)
  set(${PROJECT_NAME}.BUILD_DOCS ${DOXYGEN_FOUND})
endif()

if(${PROJECT_NAME}.MPI)
  include (loadMPI)
endif()

if(${PROJECT_NAME}.OpenMP)
  include (loadOpenMP)
endif()

if(${PROJECT_NAME}.OpenACC)
  include (loadOpenACC)
endif()

if(${PROJECT_NAME}.PACKAGE)
  include (CPack-GSMConfig)
endif()

if(${PROJECT_NAME}.SIM_DIR)
  include (SIM_DIR)
endif()

if(${PROJECT_NAME}.TEST)
  include (CTest-GSMConfig)
endif()

# ============================================================================ #

