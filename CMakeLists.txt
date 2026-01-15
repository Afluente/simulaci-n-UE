cmake_minimum_required(VERSION 3.14)
project(SimulacionUE VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Header-only library
add_library(ue_framework INTERFACE)
target_include_directories(ue_framework INTERFACE ${CMAKE_SOURCE_DIR}/src)

# Main executable
add_executable(simulacion_ue main.cpp)
target_link_libraries(simulacion_ue PRIVATE ue_framework)

# Enable warnings
target_compile_options(simulacion_ue PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>
)
