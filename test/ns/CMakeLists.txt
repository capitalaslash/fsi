include(FSI_copy_file)

FSI_copy_file(param_file_axisym param.dat)

add_executable(axisym_test axisym_test.cpp)
target_link_libraries(axisym_test ${LIBS})
