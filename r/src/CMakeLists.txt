cmake_minimum_required(VERSION 3.1)
project(tekkamaki
  VERSION 0.1.10
  LANGUAGES CXX)

find_package(blackthunnus ${PROJECT_VERSION} REQUIRED CONFIG)

include(GNUInstallDirs)
get_target_property(CXX_INCLUDEDIR blackthunnus::blackthunnus INTERFACE_INCLUDE_DIRECTORIES)
message(STATUS "INCLUDEDIR: ${CXX_INCLUDEDIR}")
get_filename_component(CXX_LIBDIR ${CXX_INCLUDEDIR}/../${CMAKE_INSTALL_LIBDIR} ABSOLUTE)
message(STATUS "LIBDIR: ${CXX_LIBDIR}")

message(STATUS "Configure Makevars")
set(CPPFLAGS "-I${CXX_INCLUDEDIR}")
set(LDFLAGS "-L${CXX_LIBDIR}")
set(LDLIBS "-lblackthunnus")
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Makevars.in
  ${CMAKE_CURRENT_SOURCE_DIR}/Makevars @ONLY
)
