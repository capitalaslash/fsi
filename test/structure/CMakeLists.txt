
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
