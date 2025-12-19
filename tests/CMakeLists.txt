# tests/CMakeLists.txt

include_directories(
    ${CMAKE_SOURCE_DIR}/src
)

set(TEST_SOURCES
    ${CMAKE_SOURCE_DIR}/src/core/Generator.cpp
    ${CMAKE_SOURCE_DIR}/src/core/SullivanFlux.cpp
    ${CMAKE_SOURCE_DIR}/src/io/RootWriter.cpp
    ${CMAKE_SOURCE_DIR}/src/models/ToyKaonModel.cpp
)

add_executable(test_run
    test_run.cpp
    ${TEST_SOURCES}
)

target_link_libraries(test_run ${ROOT_LIBRARIES})
