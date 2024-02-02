execute_process(
  COMMAND git describe --always --dirty
  OUTPUT_VARIABLE GIT_DESCRIBE
  ERROR_QUIET
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(GIT_DESCRIBE)
  set(PROJECT_VERSION ${GIT_DESCRIBE})
endif()
cmake_print_variables(PROJECT_VERSION)

file(READ "${PROJECT_SOURCE_DIR}/util/parameters.json" PARAMETERS_JSON)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/config.cpp.in
  ${CMAKE_CURRENT_SOURCE_DIR}/config.cpp @ONLY
)

# Be patient until 3.13 is popularized
target_sources(${PROJECT_NAME} PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/config.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/individual.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/population.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/program.cpp
)
