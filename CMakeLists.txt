cmake_minimum_required(VERSION 3.12)
project(kaon_sullivan_gen LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find ROOT
find_package(ROOT REQUIRED COMPONENTS RIO Tree MathCore Hist)
if (ROOT_FOUND)
  message(STATUS "Found ROOT: ${ROOT_VERSION}")
  include(${ROOT_USE_FILE})
else()
  message(FATAL_ERROR "ROOT not found")
endif()

# YAML-CPP
find_package(yaml-cpp REQUIRED)
if (NOT yaml-cpp_FOUND)
  message(FATAL_ERROR "yaml-cpp not found")
endif()

# Add subdirectories
add_subdirectory(src)
