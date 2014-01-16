
if("${CMAKE_BUILD_TYPE}" STREQUAL "")
  set(CMAKE_BUILD_TYPE NONE)
endif()
string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE_UPPER)
if(${BUILD_TYPE_UPPER} MATCHES DEBUG)
  set(METHOD dbg)
else()
  set(METHOD opt)
endif()
message(STATUS "linking against ${METHOD} libmesh library")

find_package(PkgConfig QUIET)

set(LIBMESH_DIR NOTFOUND CACHE PATH "Libmesh installtion directory")

if(LIBMESH_DIR)
  set(ENV{PKG_CONFIG_PATH} ${LIBMESH_DIR}/etc/libmesh)
endif()

pkg_check_modules(PC_LIBMESH QUIET libmesh-${METHOD})
#message(STATUS "PC_LIBMESH_FOUND = ${PC_LIBMESH_FOUND}")
#message(STATUS "PC_LIBMESH_LIBRARIES = ${PC_LIBMESH_LIBRARIES}")
#message(STATUS "PC_LIBMESH_LIBRARY_DIRS = ${PC_LIBMESH_LIBRARY_DIRS}")

# strip flags that are not definitions (-D...)
#message(STATUS "PC_LIBMESH_CFLAGS_OTHER = ${PC_LIBMESH_CFLAGS_OTHER}")
foreach(FLAG ${PC_LIBMESH_CFLAGS_OTHER})
  #message(${FLAG})
  if(${FLAG} MATCHES "^[-][D].+")
    list(APPEND ${PC_LIBMESH_CFLAGS_STRIPPED} ${FLAG})
  endif()
endforeach()
#message(STATUS "PC_LIBMESH_CFLAGS_STRIPPED = ${PC_LIBMESH_CFLAGS_STRIPPED}")
set(LIBMESH_DEFINITIONS ${PC_LIBMESH_CFLAGS_STRIPPED})

find_path(LIBMESH_INCLUDE_DIR libmesh/libmesh.h
  HINTS ${PC_LIBMESH_INCLUDEDIR} ${PC_LIBMESH_INCLUDE_DIRS}
  PATH_SUFFIXES libmesh
)

find_library(LIBMESH_LIBRARY
  NAMES mesh_${METHOD} libmesh
  HINTS ${PC_LIBMESH_LIBDIR} ${PC_LIBMESH_LIBARY_DIRS}
)

set(LIBMESH_LIBRARIES ${LIBMESH_LIBRARY})
set(LIBMESH_INCLUDE_DIRS ${LIBMESH_INCLUDE_DIR})

find_program( LIBMESH_CONFIG_EXECUTABLE
    NAMES libmesh-config
    HINTS ${LIBMESH_DIR}
    PATH_SUFFIXES bin
    DOC "libmesh-config executable" )
mark_as_advanced( LIBMESH_CONFIG_EXECUTABLE )

exec_program( ${LIBMESH_CONFIG_EXECUTABLE}
  ARGS --include
  OUTPUT_VARIABLE LMC_INC_FLAG
  RETURN_VALUE LMC_INC_RET
)
string(REPLACE " " ";" LMC_INC_LIST ${LMC_INC_FLAG})
foreach( IPATH ${LMC_INC_LIST} )
  string(REGEX REPLACE "^-I" "" IPATH ${IPATH})
  string(REGEX REPLACE "//" "/" IPATH ${IPATH})
  list(APPEND LM_INC ${IPATH})
endforeach()
set(LIBMESH_INCLUDE_DIRS ${LM_INC})


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LIBMESH DEFAULT_MSG
  LIBMESH_LIBRARY LIBMESH_INCLUDE_DIR
)

mark_as_advanced(LIBMESH_INCLUDE_DIR LIBMESH_LIBRARY)

