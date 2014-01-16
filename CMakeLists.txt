cmake_minimum_required(VERSION 2.8)

project(FSI CXX)

### cmake settings ###

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wpedantic -Wdisabled-optimization ${CMAKE_CXX_FLAGS}")

### mpi ###

find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
#add_definitions(${MPI_CXX_COMPILE_FLAGS})
#include_directories(${MPI_CXX_INCLUDE_DIRS})
#link_directories(${MPI_CXX_LIBRARY_DIRS})

### libmesh ###

find_package(LIBMESH REQUIRED)
add_definitions(${LIBMESH_DEFINITIONS})

#string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE_UPPER)
#if(BUILD_TYPE_UPPER MATCHES DEBUG)
#  add_definitions(-D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC)
#endif()

#message(STATUS "LIBMESH_INCLUDE_DIRS = ${LIBMESH_INCLUDE_DIRS}")
include_directories(${LIBMESH_INCLUDE_DIRS})

#message(STATUS "LIBMESH_LIBRARIES = ${LIBMESH_LIBRARIES}")
set(LIBS ${LIBS} ${LIBMESH_LIBRARIES})

### libxml2 (optional) ###

find_package(LibXml2)

if(LIBXML2_FOUND)
  option(FSI_HAS_LIBXML2 "use libxml2 library" ON)
else()
  option(FSI_HAS_LIBXML2 "use libxml2 library" OFF)
endif()

if(FSI_HAS_LIBXML2)
  add_definitions(${LIBXML2_DEFINITIONS})
  include_directories(${LIBXML2_INCLUDE_DIR})
  set(LIBS ${LIBS} ${LIBXML2_LIBRARIES})
endif()

### config header ###

configure_file(
  "${CMAKE_SOURCE_DIR}/FSIconfig.h.in"
  "${CMAKE_BINARY_DIR}/FSIconfig.h"
)
install(FILES ${CMAKE_BINARY_DIR}/FSIconfig.h
  DESTINATION include
)
include_directories(${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR})

### fsibase ###

set(HDRS
  laplacian.hpp
)

set(SRCS
  laplacian.cpp
)

add_library(fsibase ${SRCS})
install(TARGETS fsibase
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)
set(LIBS ${LIBS} fsibase)

### fsiassembly ###

add_subdirectory(assembly)
set(LIBS ${LIBS} fsiassembly)

### fsiutil ###

add_subdirectory(util)
set(LIBS ${LIBS} fsiutil)

### print config ###

message(STATUS "final set of libraries: ${LIBS}")
message(STATUS "CXXFLAGS:               ${CMAKE_CXX_FLAGS}")

add_subdirectory(test)

