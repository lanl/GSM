# ============================================================================ #
#
# CMakeLists.txt for the Modified_DCM sub-model library
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)
PreventInSourceBuild()

# Set sub-model's name
set(model_name "Modified_DCM")

# ============================================================================ #

# Load the sources used for the sub-model
set(sources
  ${CMAKE_CURRENT_LIST_DIR}/Parameters/Parameters.f90
  ${CMAKE_CURRENT_LIST_DIR}/Data/PhysicsData.f90
  ${CMAKE_CURRENT_LIST_DIR}/ModelData/modifiedDCMDataClass.f90
  ${CMAKE_CURRENT_LIST_DIR}/ChildTypes/mDCMChildTypes.f90
  ${CMAKE_CURRENT_LIST_DIR}/modifiedDCMClass.f90
  ${CMAKE_CURRENT_LIST_DIR}/flavor.f90
  ${CMAKE_CURRENT_LIST_DIR}/LAQ2.f
  ${CMAKE_CURRENT_LIST_DIR}/tallyFragment.f90
  ${CMAKE_CURRENT_LIST_DIR}/cascaw.f90
  ${CMAKE_CURRENT_LIST_DIR}/LAQ1.f90
  ${CMAKE_CURRENT_LIST_DIR}/upacow.f90
  ${CMAKE_CURRENT_LIST_DIR}/pauli.f90
  ${CMAKE_CURRENT_LIST_DIR}/cenum1.f90
  ${CMAKE_CURRENT_LIST_DIR}/rapid.f90
  ${CMAKE_CURRENT_LIST_DIR}/collectFragment.f90
  ${CMAKE_CURRENT_LIST_DIR}/RNDM.f90
  ${CMAKE_CURRENT_LIST_DIR}/mDCMData1.f90
  )

# List the sub-model's dependencies
set(dependencies
  General
  Random_Number_Generator
  Contracts
  Standard_DCM
  )

# ============================================================================ #

# Create the sub-model's library
add_sublibrary(${model_name} "${sources}" "${dependencies}")

# ============================================================================ #
