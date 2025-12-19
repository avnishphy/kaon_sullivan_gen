# src/CMakeLists.txt
# Build the main generator executable

# Include project source dirs
include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/src/core)
include_directories(${CMAKE_SOURCE_DIR}/src/io)
include_directories(${CMAKE_SOURCE_DIR}/src/models)
include_directories(${CMAKE_SOURCE_DIR}/src/physics)
include_directories(${CMAKE_SOURCE_DIR}/src/adapters)

set(SRC_MAIN
  main.cpp
  core/Generator.cpp
  core/SullivanFlux.cpp
  io/RootWriter.cpp
#   models/ToyKaonModel.cpp
  models/LHAPDFPdfModel.cpp
  physics/StructureFunctions.cpp
  physics/CrossSection.cpp
  adapters/APFELAdapter.cpp
  adapters/PARTONSAdapter.cpp
)

add_executable(kaon_sullivan_gen ${SRC_MAIN})

# Link libraries
target_link_libraries(kaon_sullivan_gen
  PRIVATE
    ${ROOT_LIBRARIES}
    yaml-cpp
)

# Include ROOT headers for targets
target_include_directories(kaon_sullivan_gen PRIVATE ${ROOT_INCLUDE_DIRS})

# ----------------------------
# LHAPDF (via lhapdf-config)
# ----------------------------
find_program(LHAPDF_CONFIG lhapdf-config REQUIRED)

execute_process(
    COMMAND ${LHAPDF_CONFIG} --cflags
    OUTPUT_VARIABLE LHAPDF_CFLAGS
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
    COMMAND ${LHAPDF_CONFIG} --libs
    OUTPUT_VARIABLE LHAPDF_LIBS
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

separate_arguments(LHAPDF_CFLAGS)
separate_arguments(LHAPDF_LIBS)

target_compile_options(kaon_sullivan_gen PRIVATE ${LHAPDF_CFLAGS})
target_link_libraries(kaon_sullivan_gen PRIVATE ${LHAPDF_LIBS})
