cmake_minimum_required(VERSION 3.9)
project(bench C)

add_definitions(-D_XOPEN_SOURCE=700)
add_executable(${PROJECT_NAME} ${PROJECT_SOURCE_DIR}/bench.c)
target_link_libraries(${PROJECT_NAME} delaunator-c m)
