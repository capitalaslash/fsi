
set(SRC_COMMON
  assemble_disp.cpp
  move_mesh.cpp
)

add_executable(fsi_test fsi_test.cpp ${SRC_COMMON})
target_link_libraries(fsi_test ${LIBS})

FSI_Copy_File(fsi_ibc_param param_fsi_ibc.dat)

add_executable(fsi_ibc fsi_ibc.cpp ${SRC_COMMON})
target_link_libraries(fsi_ibc ${LIBS})

