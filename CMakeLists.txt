cmake_minimum_required(VERSION 3.17)
project(TriangleCentrality-in-GraphBLAS C)

set(CMAKE_C_STANDARD 11)
include_directories(/usr/local/include)
link_directories(/usr/local/lib)
add_executable(TriangleCentrality-in-GraphBLAS src/main.c include/utils.h src/triangleCentrality.c)
target_link_libraries(TriangleCentrality-in-GraphBLAS graphblas)