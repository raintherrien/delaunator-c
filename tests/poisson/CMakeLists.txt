cmake_minimum_required(VERSION 3.9)
project(poisson C)

add_executable(${PROJECT_NAME} ${PROJECT_SOURCE_DIR}/poisson.c)
target_link_libraries(${PROJECT_NAME} delaunator-c m)
