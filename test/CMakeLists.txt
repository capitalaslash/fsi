add_executable(basic_test basic_test.cpp)
target_link_libraries(basic_test ${LIBMESH_LIBRARIES})

#add_executable(stokes systems_of_equations_ex1.C)
#target_link_libraries(stokes ${LIBMESH_LIBRARIES})

add_executable(monomial_2d monomial_2d.cpp)
target_link_libraries(monomial_2d ${LIBMESH_LIBRARIES})

add_executable(bidomain bidomain.cpp)
target_link_libraries(bidomain ${LIBMESH_LIBRARIES})

add_executable(fsi_test fsi_test.cpp)
target_link_libraries(fsi_test ${LIBS})

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
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/param.dat
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/param.dat
  COMMAND ${CMAKE_COMMAND} ARGS -E copy ${CMAKE_CURRENT_SOURCE_DIR}/param.dat ${CMAKE_CURRENT_BINARY_DIR}/param.dat
)

add_custom_target( param_file ALL
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/param.dat
)
