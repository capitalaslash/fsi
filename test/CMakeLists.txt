include(FSI_copy_file)

FSI_copy_file(param_file param.dat)
FSI_copy_file(mesh_quad-2 mesh/quad-2.msh)
FSI_copy_file(mesh_quad-22 mesh/quad-22.msh)
FSI_copy_file(mesh_quad-2x10 mesh/quad-2x10.msh)
FSI_copy_file(mesh_quad-2splx10 mesh/quad-2splx10.msh)

add_subdirectory(basic)
add_subdirectory(structure)
add_subdirectory(ns)
add_subdirectory(fsi)
