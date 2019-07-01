cmake_minimum_required(VERSION 2.8.0)

project(vgprob)

set(CXXOPTS_INCLUDE ${PROJECT_SOURCE_DIR}/deps/cxxopts/include)
set(GBWT_INCLUDE ${PROJECT_SOURCE_DIR}/deps/gbwt/include)
set(GSSW_INCLUDE ${PROJECT_SOURCE_DIR}/deps/gssw/src)
set(LIBVGIO_INCLUDE ${PROJECT_SOURCE_DIR}/deps/libvgio/include)
set(SDSL_INCLUDE ${PROJECT_SOURCE_DIR}/deps/sdsl-lite/include)

set(GBWT_LIB ${PROJECT_SOURCE_DIR}/deps/gbwt)
set(GSSW_LIB ${PROJECT_SOURCE_DIR}/deps/gssw/lib)
set(LIBVGIO_LIB ${PROJECT_SOURCE_DIR}/deps/libvgio/build)
set(SDSL_LIB ${PROJECT_SOURCE_DIR}/deps/sdsl-lite/lib)

include_directories(${CXXOPTS_INCLUDE} ${GBWT_INCLUDE} ${GSSW_INCLUDE} ${LIBVGIO_INCLUDE} ${SDSL_INCLUDE})
link_directories(${GBWT_LIB} ${GSSW_LIB} ${LIBVGIO_LIB} ${SDSL_LIB})

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

	SET(CMAKE_CXX_FLAGS "--std=c++14 -Xpreprocessor -fopenmp -g -O3")

elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")

	SET(CMAKE_CXX_FLAGS "--std=c++14 -fopenmp -g -O3")

endif()

add_executable(${PROJECT_NAME} src/main.cpp src/alignment_path.cpp src/io/register_libvg_io.cpp src/io/register_loader_saver_gbwt.cpp)
target_link_libraries(${PROJECT_NAME} gbwt gssw vgio sdsl omp protobuf)
