# src/CMakeLists.txt

# Include directories
include_directories(
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/core
    ${CMAKE_SOURCE_DIR}/src/models
    ${CMAKE_SOURCE_DIR}/src/io
)

set(CORE_SOURCES
    core/Generator.cpp
    core/SullivanFlux.cpp
)

set(IO_SOURCES
    io/RootWriter.cpp
)

set(MODEL_SOURCES
    models/ToyKaonModel.cpp
)

# Main executable
add_executable(kaon_sullivan_gen
    main.cpp
    ${CORE_SOURCES}
    ${IO_SOURCES}
    ${MODEL_SOURCES}
)

target_link_libraries(kaon_sullivan_gen ${ROOT_LIBRARIES} yaml-cpp)
