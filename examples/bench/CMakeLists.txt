cmake_minimum_required(VERSION 3.9)
project(bench VERSION 1 LANGUAGES C)

# Requires POSIX for nanosleep
if(UNIX)
    add_executable(bench ${PROJECT_SOURCE_DIR}/bench.c)
    target_compile_definitions(bench PRIVATE _XOPEN_SOURCE=700) # nanosleep
    target_link_libraries(bench delaunator-c m)
endif()

