#include "interface.hpp"

#include <libmesh/equation_systems.h>
#include <libmesh/explicit_system.h>
#include <libmesh/dof_map.h>
#include <libmesh/fe.h>
#include <libmesh/quadrature_gauss.h>
#include <libmesh/numeric_vector.h>

using namespace libMesh;

void InterfaceSystem::assemble()
{
    EquationSystems const & es = get_equation_systems();

    const MeshBase& mesh = get_mesh();

    const DofMap& dof_map = this->get_dof_map();

    std::vector<dof_id_type> dof_indices_var;

    DenseVector<Real> elem_value;

    subdomain_id_type const flag_s = es.parameters.get<subdomain_id_type>("flag_s");
    subdomain_id_type const flag_f = es.parameters.get<subdomain_id_type>("flag_f");

    AutoPtr<NumericVector<Real> > on_solid (NumericVector<Real>::build(es.comm()));
    AutoPtr<NumericVector<Real> > on_fluid (NumericVector<Real>::build(es.comm()));
    on_solid->init( this->n_dofs(), this->n_local_dofs(), false, PARALLEL );
    on_fluid->init( this->n_dofs(), this->n_local_dofs(), false, PARALLEL );

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    for ( ; el != end_el; ++el)
    {
        const Elem* elem = *el;

        dof_map.dof_indices (elem, dof_indices_var, 0);

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
        dof_map.dof_indices (elem, dof_indices_var, 0);

        elem_value.resize(dof_indices_var.size());
        for( uint k=0; k < dof_indices_var.size(); k++)
        {
            dof_id_type const current_dof = dof_indices_var[k];
            if((this->solution->first_local_index() <= current_dof) &&
                    (current_dof < this->solution->last_local_index()))
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
        this->solution->add_vector(elem_value, dof_indices_var);
    }

    this->solution->close();
    this->update();

    M_isAssembled = true;

}

void InterfaceSystem::mark_facets(subdomain_id_type const /*flag*/)
{
    libmesh_error();
//        const MeshBase& mesh = get_mesh();
//        uint const dim = mesh.mesh_dimension();

//        FEType fe_type = variable_type(0);

//        AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_type));
//        QGauss qface(dim-1, fe_type.default_quadrature_order());

//        MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
//        const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

//        for ( ; el != end_el; ++el)
//        {
//            const Elem* elem = *el;
//            for (uint s=0; s<elem->n_sides(); s++)
//            {
//                AutoPtr<Elem> side (elem->build_side(s));

//                const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
//                const std::vector<Real>& JxWface = fe_face->get_JxW();

//                fe_face->reinit(elem, s);

//                for (uint qp=0; qp<qface.n_points(); qp++)
//                {
//                    for (uint i=0; i<phi_face.size(); i++)
//                    {
//                    }
//                }

//            }
//        }
}

Real InterfaceSystem::length() const
{
    Real length = 0;
    const EquationSystems& es = get_equation_systems();
    const MeshBase& mesh = get_mesh();
    boundary_id_type const flag_int = es.parameters.get<boundary_id_type>("flag_int");
    subdomain_id_type const flag_s = es.parameters.get<subdomain_id_type>("flag_s");

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    for ( ; el != end_el; ++el)
    {
        Elem const* elem = *el;
        if(elem->subdomain_id() == flag_s)
        {
            for (uint s=0; s<elem->n_sides(); s++)
            {
                AutoPtr<Elem const> side (elem->build_side(s));
                if(mesh.boundary_info->has_boundary_id(elem, s, flag_int))
                {
                    length += side->length(0,1);
                }
            }
        }
    }

    comm().sum(length);

    return length;
}
