#ifndef INTERFACE_HPP
#define INTERFACE_HPP

#include "FSI.hpp"

#include <libmesh/explicit_system.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/dof_map.h>

class InterfaceSystem: public libMesh::ExplicitSystem
{
public:

    InterfaceSystem (libMesh::EquationSystems& es,
                     const std::string& name,
                     const unsigned int number):
        libMesh::ExplicitSystem(es, name, number),
        M_isAssembled(false),
        M_facetsMarked(false)
    {
        this->add_variable("interface", libMesh::SECOND, libMesh::LAGRANGE);
    }

    // build interface position in solution vector
    void assemble();

    // mark facets on interface with given flag
    void mark_facets(libMesh::subdomain_id_type const flag);
    
    // compute total length of the interface
    Real length() const;

    // set a bool vector to know if dofs on element are on interface
    void on_interface( libMesh::Elem const* e, std::vector<bool>& on_int) const;

private:
    bool M_isAssembled;
    bool M_facetsMarked;
    // boundary_id_type const M_fluidFlag;
    // boundary_id_type const M_solidFlag;
};

inline void InterfaceSystem::on_interface( libMesh::Elem const* e, std::vector<bool>& on_int) const
{
    std::vector<dof_id_type> dof_indices;
    get_dof_map().dof_indices(e, dof_indices, 0);

    on_int.resize(dof_indices.size());

    for(uint i = 0; i < dof_indices.size(); i++)
    {
        if((solution->first_local_index() <= dof_indices[i]) &&
                (dof_indices[i] < solution->last_local_index()) &&
                ((*solution)(dof_indices[i]) > 0.5) )
            on_int[i] = true;
        else
            on_int[i] = false;
    }
}

#endif // INTERFACE_HPP
