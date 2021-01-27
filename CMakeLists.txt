cmake_minimum_required(VERSION 3.13)
project(delaunator-c VERSION 1 LANGUAGES C)

set(CMAKE_C_STANDARD          99)
set(CMAKE_C_STANDARD_REQUIRED OFF)
set(CMAKE_C_EXTENSIONS        OFF)

# Try to compile release
set(DEFAULT_BUILD_TYPE "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

set(DELC_SOURCES
  ${PROJECT_SOURCE_DIR}/src/delaunator.c)
add_library(delaunator-c ${DELC_SOURCES})
target_include_directories(delaunator-c PUBLIC ${PROJECT_SOURCE_DIR}/include)

option(DELC_BUILD_SAN "Build with sanitizers" OFF)
if(DELC_BUILD_SAN)
    target_compile_options(delaunator-c PUBLIC -fsanitize=address,leak,undefined,pointer-compare,pointer-subtract -fstack-protector)
    target_link_options   (delaunator-c PUBLIC -fsanitize=address,leak,undefined,pointer-compare,pointer-subtract -fstack-protector)
endif()

option(DELC_BUILD_EXAMPLES "Build example programs" OFF)
if(DELC_BUILD_EXAMPLES)
  add_subdirectory(${PROJECT_SOURCE_DIR}/examples/bench)
  add_subdirectory(${PROJECT_SOURCE_DIR}/examples/poisson)
endif()

