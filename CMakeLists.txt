cmake_minimum_required(VERSION 3.10)

include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

project(rpvg)

find_package(OpenMP REQUIRED)

ExternalProject_Add(sdsl-lite-proj
  SOURCE_DIR "${PROJECT_SOURCE_DIR}/deps/sdsl-lite"
  BUILD_IN_SOURCE TRUE
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ./install.sh ${PROJECT_SOURCE_DIR}/deps/sdsl-lite
)

ExternalProject_Add(gbwt-proj
  SOURCE_DIR "${PROJECT_SOURCE_DIR}/deps/gbwt"
  BUILD_IN_SOURCE TRUE
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND ""
  INSTALL_COMMAND ""
  DEPENDS sdsl-lite-proj
)

ExternalProject_Add(gssw-proj
  SOURCE_DIR "${PROJECT_SOURCE_DIR}/deps/gssw"
  BUILD_IN_SOURCE TRUE
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND ""
  INSTALL_COMMAND ""
)

ExternalProject_Add(libvgio-proj
  SOURCE_DIR "${PROJECT_SOURCE_DIR}/deps/libvgio"
  BUILD_IN_SOURCE TRUE
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND cmake -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/deps/libvgio .
  INSTALL_COMMAND ""
)

ExternalProject_Add(libhandlegraph-proj
  SOURCE_DIR "${PROJECT_SOURCE_DIR}/deps/libhandlegraph"
  BUILD_IN_SOURCE TRUE
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND cmake -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/deps/libhandlegraph .
  INSTALL_COMMAND ""
)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

add_library(xg 
  deps/xg/src/xg.cpp
)

add_library(${PROJECT_NAME} 
  src/fragment_length_dist.cpp 
  src/paths_index.cpp
  src/alignment_path.cpp 
  src/alignment_path_finder.cpp 
  src/path_clusters.cpp 
  src/read_path_probabilities.cpp 
  src/probability_matrix_writer.cpp 
  src/path_abundance_estimator.cpp
  src/path_abundance_writer.cpp
  src/discrete_sampler.cpp
  src/io/register_libvg_io.cpp 
  src/io/register_loader_saver_gbwt.cpp 
  src/io/register_loader_saver_xg.cpp
)

include_directories(
  ${PROJECT_SOURCE_DIR}/deps/Catch2/single_include/catch2 
  ${PROJECT_SOURCE_DIR}/deps/cxxopts/include 
  ${PROJECT_SOURCE_DIR}/deps/sparsepp
  ${PROJECT_SOURCE_DIR}/deps/sdsl-lite/include
  ${PROJECT_SOURCE_DIR}/deps/gbwt/include 
  ${PROJECT_SOURCE_DIR}/deps/gssw/src 
  ${PROJECT_SOURCE_DIR}/deps/libvgio
  ${PROJECT_SOURCE_DIR}/deps/libvgio/include
  ${PROJECT_SOURCE_DIR}/deps/libhandlegraph/src/include
  ${PROJECT_SOURCE_DIR}/deps/eigen
  ${PROJECT_SOURCE_DIR}/deps/xg/src
  ${PROJECT_SOURCE_DIR}/deps/xg/deps/mmmultimap/src
  ${PROJECT_SOURCE_DIR}/deps/xg/deps/ips4o
  ${PROJECT_SOURCE_DIR}/deps/xg/deps/gfakluge/src
  ${PROJECT_SOURCE_DIR}/deps/xg/deps/gfakluge/src/tinyFA
)

link_directories(
  ${PROJECT_SOURCE_DIR}/deps/gbwt 
  ${PROJECT_SOURCE_DIR}/deps/gssw/lib 
  ${PROJECT_SOURCE_DIR}/deps/libvgio
  ${PROJECT_SOURCE_DIR}/deps/sdsl-lite/lib
  ${PROJECT_SOURCE_DIR}/deps/libhandlegraph
)

target_link_libraries(xg
  sdsl
  handlegraph
  divsufsort
  divsufsort64
)

target_link_libraries(${PROJECT_NAME}
  protobuf 
  gbwt 
  gssw 
  vgio
  sdsl
  handlegraph
  xg
  divsufsort
  divsufsort64
)

add_dependencies(xg sdsl-lite-proj)
add_dependencies(xg libhandlegraph-proj)

add_dependencies(${PROJECT_NAME} sdsl-lite-proj)
add_dependencies(${PROJECT_NAME} gbwt-proj)
add_dependencies(${PROJECT_NAME} gssw-proj)
add_dependencies(${PROJECT_NAME} libvgio-proj)
add_dependencies(${PROJECT_NAME} libhandlegraph-proj)
add_dependencies(${PROJECT_NAME} xg)

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

  set(CMAKE_CXX_FLAGS "--std=c++14 -Xpreprocessor -fopenmp -g -O3")
  target_link_libraries(${PROJECT_NAME} omp)
  target_link_libraries(xg omp)

elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")

  set(CMAKE_CXX_FLAGS "--std=c++14 -fopenmp -g -O3")
  target_link_libraries(${PROJECT_NAME} atomic)
  target_link_libraries(xg atomic)

endif()

add_executable(${PROJECT_NAME}-bin src/main.cpp) 
target_link_libraries(${PROJECT_NAME}-bin ${PROJECT_NAME})

set_target_properties(${PROJECT_NAME}-bin PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
 
add_executable(${PROJECT_NAME}-tests
  src/tests/main_test.cpp 
  src/tests/fragment_length_dist_test.cpp
  src/tests/paths_index_test.cpp
  src/tests/alignment_path_test.cpp
  src/tests/alignment_path_finder_test.cpp
  src/tests/read_path_probabilities_test.cpp
  src/tests/path_clusters_test.cpp
  src/tests/path_abundance_estimator_test.cpp
)

include_directories(
  ${PROJECT_SOURCE_DIR}/deps/catch2/single_include
)
 
target_link_libraries(${PROJECT_NAME}-tests
  ${PROJECT_NAME}
)

enable_testing()
 
add_test(
  NAME ${PROJECT_NAME}-tests
  COMMAND ${PROJECT_NAME}-tests
)

