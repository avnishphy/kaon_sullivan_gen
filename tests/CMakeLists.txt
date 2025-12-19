# tests/CMakeLists.txt
add_executable(test_run
  test_run.cpp
  ../src/core/Generator.cpp
  ../src/core/SullivanFlux.cpp
  ../src/io/RootWriter.cpp
#   ../src/models/ToyKaonModel.cpp
)

target_link_libraries(test_run PRIVATE ${ROOT_LIBRARIES} yaml-cpp)
target_include_directories(test_run PRIVATE ${CMAKE_SOURCE_DIR}/src)
