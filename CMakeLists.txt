cmake_minimum_required(VERSION 3.1)
project(blackthunnus CXX)
option(PROFILE "Use gperftools" FALSE)

include_directories(~/local/include)
include_directories(/usr/local/include)
link_directories(~/local/lib)
link_directories(/usr/local/lib)

set(CMAKE_VERBOSE_MAKEFILE TRUE)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_EXTENSIONS FALSE)
set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -Wpedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG")
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

file(READ "${CMAKE_SOURCE_DIR}/util/parameters.json" PARAMETERS_JSON)
configure_file(
  ${CMAKE_SOURCE_DIR}/src/config.hpp.in
  ${CMAKE_SOURCE_DIR}/src/config.hpp @ONLY
)

set(Boost_NO_BOOST_CMAKE TRUE)
find_package(Boost REQUIRED COMPONENTS program_options iostreams filesystem system)
message("Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
message("Boost_LIBRARY_DIRS: ${Boost_LIBRARY_DIRS}")

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "~/local" CACHE PATH "install prefix" FORCE)
endif()
message("CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")

aux_source_directory(src srcfiles)
message("Source files: ${srcfiles}")

add_executable(${PROJECT_NAME} main.cpp ${srcfiles})
target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES} sfmt)

if(PROFILE)
  target_link_libraries(${PROJECT_NAME} profiler)
endif()

install(TARGETS ${PROJECT_NAME}
  RUNTIME DESTINATION bin
)
