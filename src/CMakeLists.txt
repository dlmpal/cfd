cmake_minimum_required(VERSION 3.18)
project(hclib VERSION "0.0.1")

# Build type
if (CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE RELEASE)
endif()

# Sources
add_library(hclib SHARED
${CMAKE_CURRENT_SOURCE_DIR}/grid.cpp
${CMAKE_CURRENT_SOURCE_DIR}/grid_array.cpp
${CMAKE_CURRENT_SOURCE_DIR}/flux_function.cpp
${CMAKE_CURRENT_SOURCE_DIR}/flux_assembly.cpp
${CMAKE_CURRENT_SOURCE_DIR}/numerical_flux.cpp
${CMAKE_CURRENT_SOURCE_DIR}/equation_of_state.cpp
${CMAKE_CURRENT_SOURCE_DIR}/erk.cpp
${CMAKE_CURRENT_SOURCE_DIR}/cfl_utils.cpp
${CMAKE_CURRENT_SOURCE_DIR}/vtk.cpp)

# Release flags
if(CMAKE_BUILD_TYPE MATCHES RELEASE)
    list(APPEND cxx_release_flags -O3)
endif()

# Debug flags
if (CMAKE_BUILD_TYPE MATCHES DEBUG)
    list(APPEND  cxx_debug_flags -g;-Wall;-Wextra;-Wpedantic;-Wconversion)
endif()

target_compile_options(hclib PUBLIC "-std=c++23")

# Set flags (depending on build type)
target_compile_options(hclib PRIVATE
                       $<$<CONFIG:RELEASE>:${cxx_release_flags}>
                       $<$<CONFIG:DEBUG>:${cxx_debug_flags}>)
