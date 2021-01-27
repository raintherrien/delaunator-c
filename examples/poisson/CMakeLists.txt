cmake_minimum_required(VERSION 3.9)
project(poisson VERSION 1 LANGUAGES C)

add_executable(poisson ${PROJECT_SOURCE_DIR}/poisson.c)
target_link_libraries(poisson delaunator-c m)

