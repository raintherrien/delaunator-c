cmake_minimum_required(VERSION 3.13)
project(delaunay VERSION 1 LANGUAGES C)

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

set(DELC_SOURCES ${PROJECT_SOURCE_DIR}/src/delaunay.c)
add_library(delaunay ${DELC_SOURCES})
target_include_directories(delaunay PUBLIC
                           $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
                           $<INSTALL_INTERFACE:include>)

install(TARGETS delaunay
        EXPORT delaunay-config
		LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
		ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
		RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
		PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES ${PROJECT_SOURCE_DIR}/include/delaunay/delaunay.h
              ${PROJECT_SOURCE_DIR}/include/delaunay/helper.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/include/delaunay)
install(EXPORT delaunay-config
        FILE delaunay-config.cmake
		DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/delaunay)

option(DELC_BUILD_SAN "Build with sanitizers" OFF)
if(DELC_BUILD_SAN)
	target_compile_options(delaunay PUBLIC -fsanitize=address,leak,undefined,pointer-compare,pointer-subtract -fstack-protector)
	target_link_options   (delaunay PUBLIC -fsanitize=address,leak,undefined,pointer-compare,pointer-subtract -fstack-protector)
endif()

option(DELC_BUILD_EXAMPLES "Build example programs" OFF)
if(DELC_BUILD_EXAMPLES)
	add_subdirectory(${PROJECT_SOURCE_DIR}/examples/bench)
	add_subdirectory(${PROJECT_SOURCE_DIR}/examples/poisson)
endif()

