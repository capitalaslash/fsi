cmake_minimum_required(VERSION 2.8)

project(LibMeshApp CXX)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

set(CMAKE_CXX_FLAGS " -Wall -Wpedantic -std=c++11 ${CMAKE_CXX_FLAGS}")

find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
#add_definitions(${MPI_CXX_COMPILE_FLAGS})
#include_directories(${MPI_CXX_INCLUDE_DIRS})
#link_directories(${MPI_CXX_LIBRARY_DIRS})

find_package(LIBMESH REQUIRED)
message(STATUS "libmesh_def = ${LIBMESH_DEFINITIONS}")
add_definitions(${LIBMESH_DEFINITIONS})

#string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE_UPPER)
#if(BUILD_TYPE_UPPER MATCHES DEBUG)
#  add_definitions(-D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC)
#endif()

message(STATUS "LIBMESH_INCLUDE_DIRS = ${LIBMESH_INCLUDE_DIRS}")
include_directories(${LIBMESH_INCLUDE_DIRS})

add_executable(basic_test basic_test.cpp)
target_link_libraries(basic_test ${LIBMESH_LIBRARIES})

#add_executable(stokes systems_of_equations_ex1.C)
#target_link_libraries(stokes ${LIBMESH_LIBRARIES})

add_executable(monomial_2d monomial_2d.cpp)
target_link_libraries(monomial_2d ${LIBMESH_LIBRARIES})

add_executable(structure_static_2d structure_static_2d.cpp)
target_link_libraries(structure_static_2d ${LIBMESH_LIBRARIES})

add_executable(structure_static_3d structure_static_3d.cpp)
target_link_libraries(structure_static_3d ${LIBMESH_LIBRARIES})

add_executable(structure_2d structure_2d.cpp)
target_link_libraries(structure_2d ${LIBMESH_LIBRARIES})

add_executable(structure_3d structure_3d.cpp)
target_link_libraries(structure_3d ${LIBMESH_LIBRARIES})

add_executable(structure_seg2d structure_seg2d.cpp)
target_link_libraries(structure_seg2d ${LIBMESH_LIBRARIES})

add_executable(structure_seg3d structure_seg3d.cpp)
target_link_libraries(structure_seg3d ${LIBMESH_LIBRARIES})

add_custom_command(
  OUTPUT ${CMAKE_BINARY_DIR}/param.dat
  DEPENDS ${CMAKE_SOURCE_DIR}/param.dat
  COMMAND ${CMAKE_COMMAND} ARGS -E copy ${CMAKE_SOURCE_DIR}/param.dat ${CMAKE_BINARY_DIR}/param.dat
)

add_custom_target( param_files ALL
  DEPENDS ${CMAKE_BINARY_DIR}/param.dat
)
