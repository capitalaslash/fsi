
set(SRC_COMMON
  assemble_disp.cpp
  move_mesh.cpp
)

add_executable(fsi_test fsi_test.cpp ${SRC_COMMON})
target_link_libraries(fsi_test ${LIBS})
