cmake_minimum_required(VERSION 3.1)
include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

project(vgprob)

find_package(Protobuf REQUIRED)
find_package(OpenMP REQUIRED)

ExternalProject_Add(catch2-proj
  GIT_REPOSITORY "https://github.com/catchorg/Catch2.git"
  GIT_TAG "820271bf246afd4a381ff753b0cda14aeb1d0a1c"
  SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/catch2"
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
  GIT_TAG "fefa626dec310edef0c7269d65c59bb95e22f6dd"
  SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/gbwt"
  BUILD_IN_SOURCE TRUE
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND make
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
  BUILD_COMMAND make
  INSTALL_COMMAND ""
)

ExternalProject_Add(libvgio-proj
  GIT_REPOSITORY "https://github.com/vgteam/libvgio.git"
  GIT_TAG "f6f93ddeeb97a3977ad8109f16bf031da718e40a"
  SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/libvgio"
  BUILD_IN_SOURCE TRUE
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
)

include_directories(
  ${PROJECT_SOURCE_DIR}/external/catch2/include 
  ${PROJECT_SOURCE_DIR}/external/cxxopts/include 
  ${PROJECT_SOURCE_DIR}/external/sdsl-lite/include
  ${PROJECT_SOURCE_DIR}/external/gbwt/include 
  ${PROJECT_SOURCE_DIR}/external/gssw/src 
  ${PROJECT_SOURCE_DIR}/external/libvgio 
  ${PROJECT_SOURCE_DIR}/external/libvgio/include)

link_directories(
  ${PROJECT_SOURCE_DIR}/external/gbwt 
  ${PROJECT_SOURCE_DIR}/external/gssw/lib 
  ${PROJECT_SOURCE_DIR}/external/libvgio 
  ${PROJECT_SOURCE_DIR}/external/sdsl-lite/lib)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

add_executable(
  ${PROJECT_NAME} 
  src/main.cpp 
  src/alignment_path.cpp 
  src/alignment_path_finder.cpp 
  src/path_clusters.cpp 
  src/read_path_probs.cpp 
  src/io/register_libvg_io.cpp 
  src/io/register_loader_saver_gbwt.cpp)

add_dependencies(${PROJECT_NAME} catch2-proj)
add_dependencies(${PROJECT_NAME} cxxopts-proj)
add_dependencies(${PROJECT_NAME} sdsl-lite-proj)
add_dependencies(${PROJECT_NAME} gbwt-proj)
add_dependencies(${PROJECT_NAME} gssw-proj)
add_dependencies(${PROJECT_NAME} libvgio-proj)

target_link_libraries(
  ${PROJECT_NAME} 
  protobuf 
  sdsl 
  gbwt 
  gssw 
  vgio)

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

  set(CMAKE_CXX_FLAGS "--std=c++14 -Xpreprocessor -fopenmp -g -O3")
  target_link_libraries(${PROJECT_NAME} omp)

elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")

  set(CMAKE_CXX_FLAGS "--std=c++14 -fopenmp -g -O3")

endif()


