#ifndef EXT_PRESSURE_HPP
#define EXT_PRESSURE_HPP

#include "FSIconfig.h"

#include <libmesh/libmesh.h>

struct ExtPressure
{
    Real const max;
    Real const t_in;
    Real const t_out;
    Real const sigma;

    explicit ExtPressure (libMesh::Parameters const& p):
        max(p.get<Real>("pressure_max")),
        t_in(p.get<Real>("pressure_tin")),
        t_out(p.get<Real>("pressure_tout")),
        sigma(p.get<Real>("pressure_sigma"))
    {}

    Real time_evol(Real const t)
    {
        Real const x = (t - t_in) / (t_out - t_in);

        return max * 1 / (1. + std::exp(-5.*x));
        //        if (x < 0.) return 0.;
        //        else if (x > 1.) return max;
        //        return max * x;
    }

    Real operator() (Point const& p, Real const t)
    {
        return std::exp(-p(0)*p(0)/sigma) * time_evol(t);
    }

    Real grad( Point const& p, Real const t)
    {
        return -2.*p(0)*this->operator ()(p, t) / sigma;
    }
};

void compute_ext_pressure(EquationSystems& es)
{
    MeshBase const& mesh = es.get_mesh();

    ExplicitSystem& sys = es.get_system<ExplicitSystem>("ext_pressure");

    const DofMap & dof_map = sys.get_dof_map();

    std::vector<dof_id_type> dof_indices;

    ExtPressure p0 (es.parameters);

    subdomain_id_type const flag_f = es.parameters.get<subdomain_id_type>("flag_f");

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    for ( ; el != end_el; ++el)
    {
        const Elem* elem = *el;

        if (elem->subdomain_id() == flag_f)
        {
            dof_map.dof_indices (elem, dof_indices);
            for (uint i=0; i<dof_indices.size(); i++)
                sys.solution->set(dof_indices[i], p0(elem->point(i), sys.time));
        }
    }
}


#endif // EXT_PRESSURE_HPP
