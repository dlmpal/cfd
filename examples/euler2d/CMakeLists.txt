cmake_minimum_required(VERSION 3.18)
project(euler2d VERSION "0.0.1")

add_executable(euler2d main.cpp)
target_compile_options(euler2d PRIVATE "-std=c++23")
target_include_directories(euler2d PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../src)
target_link_directories(euler2d PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../build)
target_link_libraries(euler2d PRIVATE hclib)
