cmake_minimum_required(VERSION 3.12)

include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

project(rpvg)

find_package(OpenMP REQUIRED)

ExternalProject_Add(catch2-proj
  GIT_REPOSITORY "https://github.com/catchorg/Catch2.git"
  GIT_TAG "820271bf246afd4a381ff753b0cda14aeb1d0a1c"
  SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/catch2"
  BUILD_IN_SOURCE TRUE
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
)

ExternalProject_Add(cxxopts-proj
  GIT_REPOSITORY "https://github.com/jarro2783/cxxopts.git"
  GIT_TAG "3c73d91c0b04e2b59462f0a741be8c07024c1bc0"
  SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/cxxopts"
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
)

ExternalProject_Add(sdsl-lite-proj
  GIT_REPOSITORY "https://github.com/simongog/sdsl-lite.git"
  GIT_TAG "ddb0fbbc33bb183baa616f17eb48e261ac2a3672"
  SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/sdsl-lite"
  BUILD_IN_SOURCE TRUE
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ./install.sh ${PROJECT_SOURCE_DIR}/external/sdsl-lite
)

ExternalProject_Add(gbwt-proj
  GIT_REPOSITORY "https://github.com/jltsiren/gbwt.git"
  GIT_TAG "4f5dea6a3b7b63528cceae896379a5a5cb52b481"
  SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/gbwt"
  BUILD_IN_SOURCE TRUE
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND ""
  INSTALL_COMMAND ""
  DEPENDS sdsl-lite-proj
)

ExternalProject_Add(gssw-proj
  GIT_REPOSITORY "https://github.com/vgteam/gssw.git"
  GIT_TAG "79e0599d46133dc1c1c7dce451f29808c2f7e3a9"
  SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/gssw"
  BUILD_IN_SOURCE TRUE
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND ""
  INSTALL_COMMAND ""
)

ExternalProject_Add(libvgio-proj
  GIT_REPOSITORY "https://github.com/vgteam/libvgio.git"
  GIT_TAG "5f9e7cbd873ad11c82dcfddb0488e1b769d6dd9b"
  SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/libvgio"
  BUILD_IN_SOURCE TRUE
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND cmake -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/external/libvgio .
  INSTALL_COMMAND ""
)

ExternalProject_Add(libhandlegraph-proj
  GIT_REPOSITORY "https://github.com/vgteam/libhandlegraph.git"
  GIT_TAG "47791b556bb963d7f1922548cf512749e59dab86"
  SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/libhandlegraph"
  BUILD_IN_SOURCE TRUE
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND cmake -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/external/libhandlegraph .
  INSTALL_COMMAND ""
)

ExternalProject_Add(xg-proj
  GIT_REPOSITORY "https://github.com/jonassibbesen/xg.git"
  GIT_TAG "ef394b055aaa527d938997b6b4556fef4d6d05c1"
  SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/xg"
  BUILD_IN_SOURCE TRUE
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND cmake -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/external/xg -DINSTALL_DIR=${PROJECT_SOURCE_DIR}/external/xg -DCOMPILE_HANDLE=0 -DCOMPILE_SDSL=0 -Dhandlegraph_INCLUDE=${PROJECT_SOURCE_DIR}/external/libhandlegraph/src/include/ -Dhandlegraph_LIB=${PROJECT_SOURCE_DIR}/external/libhandlegraph -Dsdsl-lite_INCLUDE=${PROJECT_SOURCE_DIR}/external/sdsl-lite/include -Dsdsl-lite-divsufsort_INCLUDE=${PROJECT_SOURCE_DIR}/external/sdsl-lite/include -Dsdsl-lite_LIB=${PROJECT_SOURCE_DIR}/external/sdsl-lite/lib -Dsdsl-lite-divsufsort_LIB=${PROJECT_SOURCE_DIR}/external/sdsl-lite/lib .
  INSTALL_COMMAND ""
  DEPENDS sdsl-lite-proj
  DEPENDS libhandlegraph-proj
)


set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

add_library(${PROJECT_NAME} 
  src/fragment_length_dist.cpp 
  src/paths_index.cpp
  src/alignment_path.cpp 
  src/alignment_path_finder.cpp 
  src/path_clusters.cpp 
  src/read_path_probs.cpp 
  src/io/register_libvg_io.cpp 
  src/io/register_loader_saver_gbwt.cpp 
  src/io/register_loader_saver_xg.cpp
)

include_directories(
  ${PROJECT_SOURCE_DIR}/external/cxxopts/include 
  ${PROJECT_SOURCE_DIR}/external/sdsl-lite/include
  ${PROJECT_SOURCE_DIR}/external/gbwt/include 
  ${PROJECT_SOURCE_DIR}/external/gssw/src 
  ${PROJECT_SOURCE_DIR}/external/libvgio
  ${PROJECT_SOURCE_DIR}/external/libvgio/include
  ${PROJECT_SOURCE_DIR}/external/libhandlegraph/src/include
  ${PROJECT_SOURCE_DIR}/external/xg/src
  ${PROJECT_SOURCE_DIR}/external/xg/mmmultimap-prefix/src/mmmultimap/src
  ${PROJECT_SOURCE_DIR}/external/xg/ips4o-prefix/src/ips4o
)

link_directories(
  ${PROJECT_SOURCE_DIR}/external/gbwt 
  ${PROJECT_SOURCE_DIR}/external/gssw/lib 
  ${PROJECT_SOURCE_DIR}/external/libvgio
  ${PROJECT_SOURCE_DIR}/external/sdsl-lite/lib
  ${PROJECT_SOURCE_DIR}/external/libhandlegraph
  ${PROJECT_SOURCE_DIR}/external/xg/lib
)

add_dependencies(${PROJECT_NAME} catch2-proj)
add_dependencies(${PROJECT_NAME} cxxopts-proj)
add_dependencies(${PROJECT_NAME} sdsl-lite-proj)
add_dependencies(${PROJECT_NAME} gbwt-proj)
add_dependencies(${PROJECT_NAME} gssw-proj)
add_dependencies(${PROJECT_NAME} libvgio-proj)
add_dependencies(${PROJECT_NAME} libhandlegraph-proj)
add_dependencies(${PROJECT_NAME} xg-proj)

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

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

  set(CMAKE_CXX_FLAGS "--std=c++14 -Xpreprocessor -fopenmp -g -O3 -Wall")
  target_link_libraries(${PROJECT_NAME} omp)

elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")

  set(CMAKE_CXX_FLAGS "--std=c++14 -fopenmp -g -O3")

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
  src/tests/read_path_probs_test.cpp
  src/tests/path_clusters_test.cpp
)

include_directories(
  ${PROJECT_SOURCE_DIR}/external/catch2/single_include
)
 
target_link_libraries(${PROJECT_NAME}-tests
  ${PROJECT_NAME}
)

enable_testing()
 
add_test(
  NAME ${PROJECT_NAME}-tests
  COMMAND ${PROJECT_NAME}-tests
)

