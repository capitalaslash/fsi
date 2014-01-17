#macro to copy a file to the binary dir
macro(FSI_copy_file target filename)
  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${filename}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${filename}
    COMMAND ${CMAKE_COMMAND} ARGS -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${filename} ${CMAKE_CURRENT_BINARY_DIR}/${filename}
  )

  add_custom_target( ${target} ALL
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${filename}
  )
endmacro()

FSI_copy_file(param_file param.dat)
FSI_copy_file(mesh_quad-2 mesh/quad-2.msh)
FSI_copy_file(mesh_quad-22 mesh/quad-22.msh)
FSI_copy_file(mesh_quad-2x10 mesh/quad-2x10.msh)
FSI_copy_file(mesh_quad-2splx10 mesh/quad-2splx10.msh)

add_executable(basic_test basic_test.cpp)
target_link_libraries(basic_test ${LIBS})
add_test(test1 basic_test)

add_executable(laplacian_test laplacian_test.cpp)
target_link_libraries(laplacian_test ${LIBS})

add_executable(monomial_2d monomial_2d.cpp)
target_link_libraries(monomial_2d ${LIBS})

add_executable(bidomain bidomain.cpp)
target_link_libraries(bidomain ${LIBS})

add_executable(fsi_test fsi_test.cpp)
target_link_libraries(fsi_test ${LIBS})

add_executable(structure_static_2d structure_static_2d.cpp)
target_link_libraries(structure_static_2d ${LIBS})

add_executable(structure_static_3d structure_static_3d.cpp)
target_link_libraries(structure_static_3d ${LIBS})

add_executable(structure_2d structure_2d.cpp)
target_link_libraries(structure_2d ${LIBS})

add_executable(structure_3d structure_3d.cpp)
target_link_libraries(structure_3d ${LIBS})

add_executable(structure_seg2d structure_seg2d.cpp)
target_link_libraries(structure_seg2d ${LIBS})

add_executable(structure_seg3d structure_seg3d.cpp)
target_link_libraries(structure_seg3d ${LIBS})
