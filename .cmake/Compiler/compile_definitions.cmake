# ============================================================================ #
#
# Sets default compile defintions
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)


# Specify macros for release / debug
string(TOUPPER ${CMAKE_BUILD_TYPE} ${CMAKE_BUILD_TYPE})
set (contract_defines)
If (${CMAKE_BUILD_TYPE} STREQUAL "RELEASE")
  list(APPEND
    contract_defines
    "GSM_RELEASE"
    "CONTRACTS_LEVEL=0"
    )
else()
  list(APPEND
    contract_defines
    "GSM_DEBUG"
    "CONTRACTS_LEVEL=4"
    )
endif()


# Create compile definitions variable
set(${PROJECT_NAME}_compile_definitions)
list(APPEND
  ${PROJECT_NAME}_compile_definitions
  ${contract_defines}
  )


# Cache compile definitions
set(
  ${PROJECT_NAME}_compile_definitions
  ${${PROJECT_NAME}_compile_definitions}
  CACHE INTERNAL
  "Global ${PROJECT_NAME} library compile definitions"
  )
