# Configure CMake project
cmake_minimum_required(VERSION 3.1)
project(apr VERSION 0.0.0.1)

# Require C++11 standard features
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add header files directory
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")

# Add source files directory
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/src")
