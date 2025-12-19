cmake_minimum_required(VERSION 3.10)

project(kaon_sullivan_gen LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find ROOT
find_package(ROOT REQUIRED COMPONENTS Tree MathCore Hist RIO yaml-cpp)
include(${ROOT_USE_FILE})

# Add subdirectories
add_subdirectory(src)
add_subdirectory(tests)
