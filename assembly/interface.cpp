#include "interface.hpp"

#include <libmesh/explicit_system.h>
#include <libmesh/dof_map.h>
#include <libmesh/fe.h>
#include <libmesh/quadrature_gauss.h>
#include <libmesh/numeric_vector.h>

using namespace libMesh;

void assemble_interface( EquationSystems & es, std::string const & name )
{
    if (!(name == "interface"))
    {
        libmesh_error();
    }

    const MeshBase& mesh = es.get_mesh();

    const uint dim = mesh.mesh_dimension();

    ExplicitSystem & sys = es.get_system<ExplicitSystem> (name);

    const uint var = sys.variable_number(name);

    const DofMap& dof_map = sys.get_dof_map();

    FEType fe_type = sys.variable_type(var);

    AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));
    QGauss qrule (dim, fe_type.default_quadrature_order());
    fe->attach_quadrature_rule (&qrule);

    std::vector<dof_id_type> dof_indices_var;

    DenseVector<Real> elem_value;

    uint const flag_s = es.parameters.get<uint>("flag_s");
    uint const flag_f = es.parameters.get<uint>("flag_f");

    AutoPtr<NumericVector<Real> > on_solid (NumericVector<Real>::build(es.comm()));
    AutoPtr<NumericVector<Real> > on_fluid (NumericVector<Real>::build(es.comm()));
    on_solid->init( sys.n_dofs(), sys.n_local_dofs(), false, PARALLEL );
    on_fluid->init( sys.n_dofs(), sys.n_local_dofs(), false, PARALLEL );

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    for ( ; el != end_el; ++el)
    {
        const Elem* elem = *el;

        fe->reinit(elem);

        dof_map.dof_indices (elem, dof_indices_var, var);

        elem_value.resize(dof_indices_var.size());
        for( uint k=0; k < dof_indices_var.size(); k++)
        {
            elem_value(k) = 1.;
        }

        if(elem->subdomain_id() == flag_s)
        {
            on_solid->add_vector(elem_value, dof_indices_var);
        }
        else if(elem->subdomain_id() == flag_f)
        {
            on_fluid->add_vector(elem_value, dof_indices_var);
        }
        else
        {
            std::abort();
        }
    }
    on_solid->close();
    on_fluid->close();

    std::map<dof_id_type,bool> is_updated;

    el = mesh.active_local_elements_begin();

    for( ; el != end_el; ++el )
    {
        const Elem* elem = *el;

        fe->reinit(elem);

        dof_map.dof_indices (elem, dof_indices_var, var);

        elem_value.resize(dof_indices_var.size());
        for( uint k=0; k < dof_indices_var.size(); k++)
        {
            dof_id_type const current_dof = dof_indices_var[k];
            if((sys.solution->first_local_index() <= current_dof) &&
                    (current_dof < sys.solution->last_local_index()))
            {
                bool const node_on_solid = (*on_solid)(current_dof) > 0.;
                bool const node_on_fluid = (*on_fluid)(current_dof) > 0.;
                if ((!is_updated[current_dof]) &&
                        node_on_solid &&
                        node_on_fluid)
                {
                    elem_value(k) = 1.;
                    is_updated[dof_indices_var[k]] = true;
                }
            }
        }
        sys.solution->add_vector(elem_value, dof_indices_var);
    }

    sys.solution->close();
    sys.update();
}
