#include "mat.hpp"

#include <libmesh/explicit_system.h>
#include <libmesh/dof_map.h>
#include <libmesh/fe.h>
#include <libmesh/quadrature_gauss.h>
#include <libmesh/numeric_vector.h>

using namespace libMesh;

void assemble_mat( EquationSystems & es, std::string const & name )
{
    if (!(name == "mat"))
    {
        libmesh_error();
    }

    const MeshBase& mesh = es.get_mesh();

    ExplicitSystem & sys = es.get_system<ExplicitSystem> (name);

    const uint var = sys.variable_number(name);

    const DofMap& dof_map = sys.get_dof_map();

    std::vector<dof_id_type> dof_indices_var;

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    for ( ; el != end_el; ++el)
    {
        const Elem* elem = *el;

        dof_map.dof_indices (elem, dof_indices_var, var);

        // We are using CONSTANT MONOMIAL basis functions, hence we only need to get
        // one dof index per variable
        dof_id_type dof_index = dof_indices_var[0];

        if( (sys.solution->first_local_index() <= dof_index) &&
                (dof_index < sys.solution->last_local_index()))
        {
            sys.solution->set(dof_index, elem->subdomain_id());
        }

    }

    sys.solution->close();
    sys.update();
}
