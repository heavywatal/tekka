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
  config.cpp.in
  config.cpp @ONLY
)

target_sources(${PROJECT_NAME} PRIVATE
  ${CMAKE_CURRENT_BINARY_DIR}/config.cpp
  individual.cpp
  population.cpp
  program.cpp
)
