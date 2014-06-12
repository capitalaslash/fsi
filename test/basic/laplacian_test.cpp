#include "FSI.hpp"

#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/equation_systems.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/dirichlet_boundaries.h>
#include <libmesh/zero_function.h>
#include <libmesh/dof_map.h>
#include <libmesh/vtk_io.h>

#include "laplacian.hpp"

int main( int argc, char** argv )
{
    LibMeshInit init (argc, argv);

    Mesh mesh(init.comm());

    MeshTools::Generation::build_square (mesh,
                                         4, 4,
                                         0.0, 1.0,
                                         0.0, 1.0,
                                         QUAD4);

//    XdrIO mesh_io(mesh);
//    mesh_io.read("one_tri.xda");

    mesh.print_info();

    EquationSystems es (mesh);

    LinearImplicitSystem& system = es.add_system<LinearImplicitSystem>("lap");

    uint u_var = system.add_variable("u", FIRST, LAGRANGE);

    Laplacian lap(es);

    system.attach_assemble_object(lap);

    std::set<boundary_id_type> bd_ids;
    bd_ids.insert(1);
    bd_ids.insert(3);

    std::vector<uint> vars(1,u_var);

    ZeroFunction<Real> zero;

    DirichletBoundary dirichlet_bc(bd_ids, vars, &zero);

    system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);

    es.init();

    es.print_info();

    system.solve();

    VTKIO(mesh).write_equation_systems("lap.pvtu",es);

    return 0;
}
