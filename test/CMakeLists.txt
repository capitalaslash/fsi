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

add_subdirectory(basic)
add_subdirectory(structure)
add_subdirectory(ns)
add_subdirectory(fsi)
