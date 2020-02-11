cmake_minimum_required(VERSION 3.12)

include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

project(rpvg)

find_package(OpenMP REQUIRED)

ExternalProject_Add(catch2-proj
  GIT_REPOSITORY "https://github.com/catchorg/Catch2.git"
  GIT_TAG "b3b07215d1ca2224aea6ff3e21d87ad0f7750df2"
  SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/catch2"
  BUILD_IN_SOURCE TRUE
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
)

ExternalProject_Add(cxxopts-proj
  GIT_REPOSITORY "https://github.com/jarro2783/cxxopts.git"
  GIT_TAG "943cdcfbbaa205f9a1789d96c7a57accfabcd233"
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
  GIT_TAG "9f8efbca6dbfa3b7a4e6a65d0ef7cf5009e4aca1"
  SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/gssw"
  BUILD_IN_SOURCE TRUE
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND ""
  INSTALL_COMMAND ""
)

ExternalProject_Add(libvgio-proj
  GIT_REPOSITORY "https://github.com/vgteam/libvgio.git"
  GIT_TAG "fc9523feb660815b02a932e8e7f58467d462455e"
  SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/libvgio"
  BUILD_IN_SOURCE TRUE
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND cmake -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/external/libvgio .
  INSTALL_COMMAND ""
)

ExternalProject_Add(libhandlegraph-proj
  GIT_REPOSITORY "https://github.com/vgteam/libhandlegraph.git"
  GIT_TAG "962018b84d8525ca93e79efd5f3726a0683bcfdf"
  SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/libhandlegraph"
  BUILD_IN_SOURCE TRUE
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND cmake -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/external/libhandlegraph .
  INSTALL_COMMAND ""
)

# ExternalProject_Add(xg-proj
#   GIT_REPOSITORY "https://github.com/jonassibbesen/xg.git"
#   GIT_TAG "aa3d221d34bf2227ba42ca0c1e189536f7b6962d"
#   SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/xg"
#   BUILD_IN_SOURCE TRUE
#   UPDATE_COMMAND ""
#   CONFIGURE_COMMAND cmake -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/external/xg -DINSTALL_DIR=${PROJECT_SOURCE_DIR}/external/xg -DIS_RPVG_LIB=1 -Dhandlegraph_INCLUDE=${PROJECT_SOURCE_DIR}/external/libhandlegraph/src/include/ -Dhandlegraph_LIB=${PROJECT_SOURCE_DIR}/external/libhandlegraph -Dsdsl-lite_INCLUDE=${PROJECT_SOURCE_DIR}/external/sdsl-lite/include -Dsdsl-lite-divsufsort_INCLUDE=${PROJECT_SOURCE_DIR}/external/sdsl-lite/include -Dsdsl-lite_LIB=${PROJECT_SOURCE_DIR}/external/sdsl-lite/lib -Dsdsl-lite-divsufsort_LIB=${PROJECT_SOURCE_DIR}/external/sdsl-lite/lib .
#   INSTALL_COMMAND ""
#   DEPENDS sdsl-lite-proj
#   DEPENDS libhandlegraph-proj
# )

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
  # src/io/register_loader_saver_xg.cpp
)

include_directories(
  ${PROJECT_SOURCE_DIR}/external/cxxopts/include 
  ${PROJECT_SOURCE_DIR}/external/sdsl-lite/include
  ${PROJECT_SOURCE_DIR}/external/gbwt/include 
  ${PROJECT_SOURCE_DIR}/external/gssw/src 
  ${PROJECT_SOURCE_DIR}/external/libvgio
  ${PROJECT_SOURCE_DIR}/external/libvgio/include
  ${PROJECT_SOURCE_DIR}/external/libhandlegraph/src/include
  # ${PROJECT_SOURCE_DIR}/external/xg/src
  # ${PROJECT_SOURCE_DIR}/external/xg/deps/mmmultimap/src
  # ${PROJECT_SOURCE_DIR}/external/xg/deps/ips4o
)

link_directories(
  ${PROJECT_SOURCE_DIR}/external/gbwt 
  ${PROJECT_SOURCE_DIR}/external/gssw/lib 
  ${PROJECT_SOURCE_DIR}/external/libvgio
  ${PROJECT_SOURCE_DIR}/external/sdsl-lite/lib
  ${PROJECT_SOURCE_DIR}/external/libhandlegraph
  # ${PROJECT_SOURCE_DIR}/external/xg/lib
)

add_dependencies(${PROJECT_NAME} catch2-proj)
add_dependencies(${PROJECT_NAME} cxxopts-proj)
add_dependencies(${PROJECT_NAME} sdsl-lite-proj)
add_dependencies(${PROJECT_NAME} gbwt-proj)
add_dependencies(${PROJECT_NAME} gssw-proj)
add_dependencies(${PROJECT_NAME} libvgio-proj)
add_dependencies(${PROJECT_NAME} libhandlegraph-proj)
# add_dependencies(${PROJECT_NAME} xg-proj)

target_link_libraries(${PROJECT_NAME}
  protobuf 
  gbwt 
  gssw 
  vgio
  sdsl
  handlegraph
  # xg
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

