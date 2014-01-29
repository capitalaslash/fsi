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
