cmake_minimum_required(VERSION 3.18)
project(euler1d VERSION "0.0.1")

add_executable(euler1d main.cpp)
target_compile_options(euler1d PRIVATE "-std=c++23")
target_include_directories(euler1d PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../src)
target_link_directories(euler1d PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../build)
target_link_libraries(euler1d PRIVATE hclib)
