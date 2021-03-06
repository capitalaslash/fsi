include(FSI_copy_file)
include(FSI_generate_mesh)

set(HDRS
  assemle_disp.hpp
  move_mesh.hpp
)

set(SRCS
  assemble_disp.cpp
  move_mesh.cpp
)

add_library(fsitest ${SRCS})

add_executable(fsi_test fsi_test.cpp)
target_link_libraries(fsi_test fsitest ${LIBS})

FSI_Copy_File(fsi_ibc_param param_fsi_ibc.dat)
add_executable(fsi_ibc fsi_ibc.cpp)
target_link_libraries(fsi_ibc fsitest ${LIBS})

FSI_Copy_File(fsi_tube_param param_fsi_tube.dat)
add_executable(fsi_tube fsi_tube.cpp)
target_link_libraries(fsi_tube fsitest ${LIBS})

FSI_Copy_File(fsi_valve_param param_fsi_valve.dat)
add_executable(fsi_valve fsi_valve.cpp)
target_link_libraries(fsi_valve fsitest ${LIBS})

FSI_generate_mesh(valve
  MESH_DIR ${CMAKE_SOURCE_DIR}/test/mesh
)
