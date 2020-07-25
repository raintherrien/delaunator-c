cmake_minimum_required(VERSION 3.9)
project(delaunator-c C)

# Options
option(DELC_BUILD_TESTS "Build test programs" OFF)
option(DELC_BUILD_SAN   "Build with sanitizers" OFF)

# C standard
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED YES)
set(CMAKE_C_EXTENSIONS NO)

# Try to compile release
set(DEFAULT_BUILD_TYPE "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Compiler flags
set(DELC_COMMON_C_FLAGS "-Wall -Werror -Wextra -pedantic")
if(DELC_BUILD_SAN)
    set(DELC_COMMON_C_FLAGS "${DELC_COMMON_C_FLAGS} -fsanitize=address -fsanitize=leak -fsanitize=undefined")
endif()
set(CMAKE_C_FLAGS_DEBUG          "${CMAKE_C_FLAGS_DEBUG}          ${DELC_COMMON_C_FLAGS} -g3")
set(CMAKE_C_FLAGS_RELEASE        "${CMAKE_C_FLAGS_RELEASE}        ${DELC_COMMON_C_FLAGS}")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} ${DELC_COMMON_C_FLAGS} -g3")
set(CMAKE_C_FLAGS_MINSIZEREL     "${CMAKE_C_FLAGS_MINSIZEREL}     ${DELC_COMMON_C_FLAGS}")

# Try to enable LTO
include(CheckIPOSupported)
check_ipo_supported(RESULT IS_IPO_SUPPORTED OUTPUT output)
if(IS_IPO_SUPPORTED)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

# Executable, source dir
set(DELC_SOURCES
  ${PROJECT_SOURCE_DIR}/src/delaunator.c)
add_library(${PROJECT_NAME} STATIC ${DELC_SOURCES})
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/include)

# Tests
if(DELC_BUILD_TESTS)
  add_subdirectory(${CMAKE_SOURCE_DIR}/tests)
endif()
