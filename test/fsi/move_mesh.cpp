
#include <libmesh/libmesh.h>
#include <libmesh/equation_systems.h>
#include <libmesh/transient_system.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/dof_map.h>

using namespace libMesh;

void move_mesh( EquationSystems& es )
{
    MeshBase& mesh = es.get_mesh();

    // MeshBase::node_iterator nd = mesh.local_nodes_begin();
    // const MeshBase::node_iterator end_nd = mesh.local_nodes_end();
    // for ( ; nd != end_nd; ++nd)
    // {
    //     Node* node = *nd;
    //     Real & x = node->operator()(0);
    //     x += .1*x;
    // }

    TransientLinearImplicitSystem & system_d =
            es.get_system<TransientLinearImplicitSystem> ("dx");
    TransientLinearImplicitSystem & system_e =
            es.get_system<TransientLinearImplicitSystem> ("dy");

    const uint d_var = system_d.variable_number ("dx");
    const uint e_var = system_e.variable_number ("dy");
    const DofMap & dof_map_d = system_d.get_dof_map();
    const DofMap & dof_map_e = system_e.get_dof_map();
    std::vector<dof_id_type> dof_indices_d;
    std::vector<dof_id_type> dof_indices_e;

    Real factor = es.parameters.get<Real> ("ale_factor");

    std::map<dof_id_type,bool> is_updated;

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    for ( ; el != end_el; ++el)
    {
        const Elem* elem = *el;

        dof_map_d.dof_indices (elem, dof_indices_d, d_var);
        dof_map_e.dof_indices (elem, dof_indices_e, e_var);

        for (uint n = 0; n < elem->n_nodes(); n++)
        {
            dof_id_type const & node_id = elem->node(n);
            if (!is_updated[node_id])
            {
                Node & node = mesh.node(node_id);
                node(0) += factor*system_d.current_solution (dof_indices_d[n]);
                node(1) += factor*system_e.current_solution (dof_indices_e[n]);
                is_updated[node_id] = true;
            }
        }
    }
}
