#macro to generate a mesh using gmsh

include(FSI_parse_opt_args)

function(FSI_generate_mesh filename)

  FSI_parse_opt_args("MESH_DIR" ${ARGN})

  if(NOT MESH_DIR)
    set(MESH_DIR ${CMAKE_CURRENT_SOURCE_DIR}/mesh)
  endif()

  add_custom_command(OUTPUT ${filename}.msh
    COMMAND gmsh ARGS -2 -order 2 ${MESH_DIR}/${filename}.geo -o ${CMAKE_CURRENT_BINARY_DIR}/${filename}.msh
  )
  add_custom_target(mesh_${filename} DEPENDS ${filename}.msh)
endfunction()

